#!/usr/bin/env python
# 2019 Mohammad Alanjary
# Wageningen University and Research (WUR)
# Bioinformatics department, Lab of Prof. Marnix Medema
#
# License: A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import os, subprocess, setlog, sqlite3 as sql
import sqlmethods
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import CompoundLocation
from Bio.Seq import Seq
from StringIO import StringIO
from modeller import *
import multiprocessing as mp

def getNRPSdomString(seqfeat):
    nrps_pks_doms = [x.split()[1] for x in seqfeat.qualifiers.get("NRPS_PKS",[]) if "Domain:" in x]
    return "->-".join(nrps_pks_doms).lower()

def hasNRPSmods(seqfeat):
    #Test that at least one C and A domain exist in the region - ASv5 annotations
    nrps_domstring = getNRPSdomString(seqfeat)
    if ('condensation' in nrps_domstring or 'cglyc' in nrps_domstring) and "amp-binding" in nrps_domstring:
        return True
    return False

def isNrps(seqfeat):
    #Looks for identifirers from antismashv4 or v5 annotations
    if any(['type: nrps' in x.lower() for x in seqfeat.qualifiers.get("sec_met",[])]):
        return True
    elif 'non-ribosomal peptide synthetase' in seqfeat.qualifiers.get("product",[]):
        return True
    elif 'peptide synthetase' in seqfeat.qualifiers.get("product",[]):
        return True
    elif 'nrps' in seqfeat.qualifiers.get("product",[]):
        return True
    elif hasNRPSmods(seqfeat):
        return True
    return False


def sortorfs(nrpsorfs):
    #Attempt to order first and last orfs based on presence of domains (TE, C-starter)
    first = []
    mid = []
    last = []
    for orf in nrpsorfs:
        #Sort based on Starter and terminator domains
        domstring = "->-".join([getNRPSdomString(x) for x in orf.features if x.type=="CDS"])
        if "condensation_starter" in domstring:
           first.append(orf)
        elif "->-td" in domstring or "thioesterase" in domstring:
            last.append(orf)
        else:
            mid.append(orf)
    return first+mid+last

def getNRPSorfs(rec):
    return [x for x in rec.features if x.type=="CDS" and isNrps(x)]

def getdomains(seq):
    #return ordered list of domains based on antismash annotations
    return [x for x in seq.features if x.type=='aSDomain']

def isCdomain(dom):
    if "Condensation" in dom.qualifiers.get('aSDomain',[""])[0] or "Cglyc" in dom.qualifiers.get('aSDomain',[""])[0] or "Heterocyclization" in dom.qualifiers.get('aSDomain',[""])[0]:
        return True
    return False

def isTerminal(dom):
    if "Thioesterase" in dom.qualifiers.get('aSDomain',[""])[0] or "TD" in dom.qualifiers.get('aSDomain',[""])[0]:
        return True
    return False

def sepmodules(domains):
    #assume ordered sequence, split modules by C-domains
    moduleList = []
    CdomainList = []
    m1start = -1
    skip = 1
    #log cdomains in each module
    if isCdomain(domains[0]):
        CdomainList.append(domains[0])
        m1start = 0
    temp = [domains[0]]
    if "NRPS-COM" in temp[0].qualifiers.get('aSDomain',[""])[0]:
        #In case NRPS-Com is before C-domain
        temp.append(domains[1])
        skip += 1
    for dom in domains[skip:]:  #skip first and add domains
        temp.append(dom)
        if isCdomain(dom):
            moduleList.append(temp)
            CdomainList.append(dom)
            temp = [dom]
    if len(temp) > 1:
        moduleList.append(temp)
    return moduleList,CdomainList,m1start

def mergelocations(seqlist,rec,tol=100):
    #merges feature locations if they are on same strand and within tolerance bp
    fuseborder = []
    #Get list of all locations and simplify complex joins with start/end
    loclist = [FeatureLocation(x.location.start,x.location.end,x.location.strand) for x in seqlist]
    domlist = [getdomains(x.extract(rec)) for x in seqlist]
    listlen = len(loclist)
    for i in range(listlen-1):
        Tdoms = [x.qualifiers.get('aSDomain',[''])[0] for x in domlist[i] if isTerminal(x)]
        if len(Tdoms):
            #Do not merge if terminal is present
            continue
        if loclist[i].strand == loclist[i+1].strand and abs(loclist[i].end-loclist[i+1].start) <= tol:
            if len(fuseborder) > 0 and fuseborder[-1][1]==i:
                #combine previous fusions
                fuseborder[-1][1]=i+1
            else:
                fuseborder.append([i,i+1])
    #revesed order to merge ending items first
    for fusion in fuseborder[::-1]:
        fused = FeatureLocation(loclist[fusion[0]].start,loclist[fusion[1]].end,loclist[fusion[0]].strand)
        #Override fused locations:
        loclist[fusion[0]:fusion[1]+1] = [fused]
    #Create list of new sequence recs for all merged orfs
    newseqs = [x.extract(rec) for x in loclist]
    #Overlapping CDS features will include domains from downstream CDS...
    # Fix this by finding all out of bounds features and removing them
    # removefeats = {}
    # for j in range(len(newseqs)):
    #     cds = [x for x in newseqs[j].features if x.type=="CDS"]
    #     if len(cds) > 1:
    #         log.debug("Overlapping CDS regions detected")
    #         for i,feat in enumerate(newseqs[j].features):
    #             if feat.location.end > cds[0].location.end or feat.location.start < cds[0].location.start:
    #                 if j not in removefeats:
    #                     removefeats[j] = []
    #                 removefeats[j].append(i)
    # for j,remlist in removefeats.items():
    #     oldlen = len(newseqs[j].features)
    #     newseqs[j].features = [f for i,f in enumerate(newseqs[j].features) if i not in remlist]
    #     log.debug("Removed %s features"%(oldlen - len(newseqs[j].features)))
    return newseqs

def mapCdomcoords(refalgn,qalign,featlocs):
    #Translates reference locations of features based on ref and query alignment
    #Featlocs: Reference locations dictionary of features relative to reference ex: {"alphahelix5":(153,165),...}
    #Get list of coordinates needing to be translated
    reflocs = [v[0] for k,v in featlocs.items() if not k.startswith("_")] + [v[1] for k,v in featlocs.items() if not k.startswith("_")]
    #Get new alignment coords for each reference coord
    j = 0 #index of ungapped position (reference position)
    newlocs = {k:None for k in reflocs} #Dictionary to translate align location
    for i,x in enumerate(refalgn):
        if x != "-": #increment ungapped position
            j += 1
        if j in reflocs: #Record alignment position in newlocs
            newlocs[j] = i+1
    #Translate the alignment locations to query locations
    L = lambda seq,i : len(seq[:i])-seq[:i].count("-")
    translocs = {k:L(qalign,v) for k,v in newlocs.items()}
    #Make new feature location dictionary
    newfeatlocs = {}
    for k,v in featlocs.items():
        if k.startswith("_"):
            continue
        newfeatlocs[k] = (translocs[v[0]],translocs[v[1]])
    return newfeatlocs

def parsePIR(alignment):
    #Read and parse the output
    if isinstance(alignment,basestring):
        with open(alignment,'r') as fil:
            result = "".join(fil.readlines())
    else: #stringIO instance
        result = alignment.getvalue()
    result = result.split(">")
    #Store alignment for Ref (first record) and query (next) based on PIR format
    refalgn = "".join(result[1].split('\n')[2:]).strip().replace("*","")
    qalgn = "".join(result[2].split('\n')[2:]).strip().replace("*","")
    return refalgn,qalgn

def parseHMM(alignments):
    #Convert stockholm
    stkhlm = AlignIO.read(StringIO(alignments),"stockholm")
    fasta = StringIO()
    status = SeqIO.write(stkhlm,fasta,"fasta")
    with open("test_hmmalign.faa","w") as fil:
        status = SeqIO.write(stkhlm,fil,"fasta")
    fasta = fasta.getvalue()
    alignList = ["".join(x.split("\n")[1:]) for x in fasta.split(">")[1:]]
    return [(alignList[0],qalign) for qalign in alignList[1:]]

def alignCdoms(cdomains,refpdb,refseq,orf,offset="FIRST",method="MDL",UseTrans=True,precompute=False,alignmentprefix=""):
    #MODELLER alignment
    alignments = []
    if method=="MDL":
        for i,Cdom in enumerate(cdomains):
            # CdomSeq = Cdom.qualifiers.get("translation",[""])[0]
            # if not CdomSeq:
            #     print "No translation!"
            newlocation = FeatureLocation(Cdom.location.start,Cdom.location.end+Cdomextend,Cdom.location.strand)
            CdomSeq = newlocation.extract(orf).seq.translate(to_stop=True)
            #format input into PIR
            inSeq = StringIO(">P1;CdomMod\nsequence:CdomMod:::::::0.00: 0.00\n%s*\n"%CdomSeq)
            if not alignmentprefix:
                output = StringIO()
            else:
                output = alignmentprefix + "_Cdom_%s.pir"%i
            #run Modeller
            env = environ()
            env.libs.topology.read('${LIB}/top_heav.lib')
            aln = alignment(env)
            templateId = os.path.splitext(os.path.split(refpdb)[-1])[0]
            mdl = model(env, file=refpdb, model_segment=('%s:A'%offset,'LAST:A'))
            aln.append_model(mdl, align_codes=templateId, atom_files=refpdb)
            aln.append(file=inSeq,align_codes="CdomMod")
            status = aln.align2d()
            aln.write(file=output, alignment_format='PIR')
            alignments.append(parsePIR(output))
    elif method=="HMM":
        #format input into FASTA
        with open(refseq,"r") as fil:
            fastaStr = "".join(fil.readlines())+"\n"
        for i,Cdom in enumerate(cdomains):
            # CdomSeq = Cdom.qualifiers.get("translation",[""])[0]
            # if not CdomSeq or UseTrans:
            #     print "Using translated sequence"

            #Use the translated seq with 150 Amino Acids added to capture latch and full C-domain
            newlocation = FeatureLocation(Cdom.location.start,Cdom.location.end+Cdomextend,Cdom.location.strand)
            CdomSeq = newlocation.extract(orf).seq.translate()
            fastaStr += ">Cdom_%s\n%s\n"%(i,CdomSeq)
        hmmFile = os.path.join(os.path.dirname(os.path.abspath(__file__)),"Condensation.hmm")
        # inSeq = "temp_%s.faa"%uuid.uuid4()
        # with open(inSeq,"w") as fil:
        #     fil.write(fastaStr)
        #Run Hmm Align
        p = subprocess.Popen(["hmmalign",hmmFile,"/dev/stdin"],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        outalign,out_err = p.communicate(input=fastaStr)
        alignments = parseHMM(outalign)
    return alignments

#simple check that A and PCP domains are only found once per module
def singledom(domlist):
    if domlist.count("PCP") == 1 and domlist.count("AMP-binding") == 1:
        return True
    elif domlist.count("AMP-binding") > 1:
        log.debug("Encountered multiple A domains")
    elif domlist.count("PCP") > 1:
        log.debug("Encountered multiple PCP domains")
    return False

def makerec(acc,clusterid,orfnum,modnum,module,Ccoords,orf,upid="",downid="",organism="",taxonomy=""):
    #Create module reccord
    modid = "%s|%s-%s"%(clusterid,orfnum,modnum)
    moddict = {"modid":modid,"modnum":modnum,"accession":acc,'clusterid':clusterid,'type':"Other","Edom":False,"domaincount":len(module),\
               "organism":organism,"taxonomy":taxonomy,"upid":upid,"downid":downid,"AdomSpec":None}
    moddict['domlist'] = [x.qualifiers.get('domain_subtype',[False])[0] if x.qualifiers.get('domain_subtype',[False])[0] else x.qualifiers.get('aSDomain',['?'])[0] for x in module]
    moddict['domstring'] = "->-".join(moddict['domlist'])
    domlocs = [x.location for x in module]
    moddict['domlocs'] = [str(x) for x in domlocs]
    locus = FeatureLocation(domlocs[0].start,domlocs[-1].end+Cdomextend,domlocs[0].strand)
    modseqrec = locus.extract(orf)
    moddict["seq"] = str(modseqrec.seq)
    moddict["domseqs"] = [str(x.extract(orf).seq) for x in domlocs]
    ##Set the type of module: Start=Start,Simple-Extender=SimpleExt,Extender=Ext,Terminal=Term,Other
    if "Starter" in moddict['domstring']:
        moddict['type'] = "CStart"
    elif moddict['domstring'] == "AMP-binding->-PCP->-Condensation" or moddict['domstring'] == "AMP-binding->-PCP->-Cglyc" :
        moddict['type'] = "Start"
    elif isCdomain(module[0]) and isCdomain(module[-1]) and singledom(moddict['domlist']):
        if len(moddict['domlist']) < 5 and moddict['domlist'][1] == "AMP-binding" and moddict['domlist'][2] == "PCP":
            moddict['type'] = "SimpleExt"
        else:
            moddict['type'] = "Ext"
    elif "->-TD" in moddict['domstring']:
        moddict['type'] = "TermD"
    elif "Thioesterase" in moddict['domstring']:
        moddict['type'] = "TermE"
    #Log tailoring and other added domains:
    moddict["extraDomains"] = [x for x in moddict['domlist'] if not x in ["Condensation","Condensation_LCL","Condensation_DCL","Condensation_Dual","Condensation_Starter","Cglyc","TD","Thioesterase","PCP","Epimerization","AMP-binding","Heterocyclization"]]
    if "Epimerization" in moddict['domstring']:
        moddict["Edom"] = True

#Store Cdomain features / sequences For C1 and C2
    def expandCdomainprops(i,C,startPos):
        for key in Ccoords[C]:
            dnalocs = tuple((x*3)+startPos for x in Ccoords[C][key]) #convert amino acid to DNA position
            moddict["%sloc#%s"%(C,key)] = list(dnalocs)
            # moddict["%sSeq#%s"%(C,key)] = moddict["seq"][moddict["%sloc#%s"%(C,key)][0]:moddict["%sloc#%s"%(C,key)][-1]]
            moddict["%sSeq#%s"%(C,key)] = moddict["seq"][moddict["%sloc#%s"%(C,key)][0]:moddict["%sloc#%s"%(C,key)][-1]]
            moddict["%sPSeq#%s"%(C,key)] = str(Seq(moddict["%sSeq#%s"%(C,key)]).translate())
        moddict["%stype"%C] = "None"
        if "_" in moddict['domlist'][i]:
            moddict["%stype"%C] = moddict['domlist'][i].split("_")[-1]
        elif "Cglyc" in moddict['domlist'][i]:
            moddict["%stype"%C] = "Cglyc"
    if isCdomain(module[0]):
        expandCdomainprops(0,"C1",0)
    if isCdomain(module[-1]):
        expandCdomainprops(-1,"C2",module[-1].location.start - module[0].location.start)
    #Store A-domain specificity
    if "AMP-binding" in moddict['domlist']:
        AdomSpec = [x.qualifiers.get('specificity',['None: None'])[0] for x in module if x.qualifiers.get('aSDomain',[''])[0]=="AMP-binding"][0]
        moddict["AdomSpecMethod"],moddict["AdomSpec"] = AdomSpec.split(": ")
    return moddict,modid

def processclustgbks_ASv5(filename,mergetol=100,alignmethod="MDL",refcoords=None,db="",savealignments=False,alignfolder="",refpdb="",refseq=""):
    results = []
    #for filename in filenames:
    #Parse the GBK file
    recs = SeqIO.parse(filename,"genbank")
    #Find all NRPS genes and extract modules
    log.info("Starting %s"%filename)
    for contignum,rec in enumerate(recs):
        clusterid = os.path.splitext(os.path.split(filename)[-1])[0]
        organism = rec.annotations.get("organism","")
        taxonomy = "|".join(rec.annotations.get("taxonomy",[""]))

        #Get NRPS orfs and merge connected, Attempt to order by start and end
        orfs = sortorfs(mergelocations(getNRPSorfs(rec),rec))
        for orfnum,orf in enumerate(orfs):
            #Get protein accession
            acc = "n/a"
            if hasattr(orf, 'id') and len(orf.id):
                acc = orf.id
            elif hasattr(orf, 'qualifiers'):
                acc = orf.qualifiers.get("protein_id",["n/a"])[0]

            alignmentprefix = ""
            if savealignments:
                alignmentprefix = os.path.join(alignfolder,"%s_%s-%s"%(clusterid,contignum,orfnum))


            orfdomains = getdomains(orf)
            if not orfdomains:
                #Skip if orf has no domains
                log.debug("No domains found, Skipping NRPS: %s|%s|%s"%(acc,clusterid,orfnum))
                continue
            modules,cdomains,m1start = sepmodules(orfdomains)
            #Skip NRPSs that dont have C domains
            if len(cdomains) < 1:
                log.debug("No C-domains found, Skipping NRPS: %s|%s|%s"%(acc,clusterid,orfnum))
                continue
            #Strore (reference align, query cdom align) for mapping
            alignments = alignCdoms(cdomains, refpdb, refseq, orf, method=alignmethod, offset=refcoords["_offset"], alignmentprefix=alignmentprefix)
            cdomCoords = [mapCdomcoords(ra,qa,refcoords) for ra,qa in alignments]

            #Make record for each module
            lastmodid = None
            uprec = None
            for i,module in enumerate(modules):
                Ccoords = {"C1":[],"C2":[]}
                #Link cdom and module indices
                c1ind = m1start + i
                c2ind = c1ind + 1
                if c1ind >= 0:
                    Ccoords["C1"] = cdomCoords[c1ind]
                if c2ind < len(cdomCoords):
                    Ccoords["C2"] = cdomCoords[c2ind]
                modrec,lastmodid = makerec(acc,clusterid,orfnum,i,module,Ccoords,orf,upid=lastmodid,organism=organism,taxonomy=taxonomy)
                #modify downstream id of last upstream record
                if uprec:
                    results[-1]["downid"] = lastmodid
                results.append(modrec)
                uprec = results[-1]
    return results
    # if db:
    #     sqlmethods.writetodb(results,db)
    # return {v["modid"]:v for v in results}

def mpASv5(filenames,mergetol=100,alignmethod="MDL",ref="2jgp",db="",savealignments=False,alignfolder="",ncpus=None):
    #Setup directories
    if not os.path.isdir(alignfolder):
        os.mkdir(alignfolder)
    #Mapped locations to reference xray structures
    refpdb = os.path.join(os.path.dirname(os.path.abspath(__file__)),"refmodels","%s.pdb"%ref)
    refseq = os.path.join(os.path.dirname(os.path.abspath(__file__)),"refmodels","%s.faa"%ref)
    coords={}
    coords['2jgp'] = {'floor': (263, 268), 'activesite': (123, 130), 'latch_loop': (338, 367), 'B1': (25,32),'acceptNlobe1': (75,123), 'acceptNlobe2': (25,49),'acceptClobe1': (332,412), 'acceptClobe2': (204,228), 'Nlobe':(0,167), 'A5': (154, 167), 'A6': (168, 183), 'Clobe':(168,423), 'floor_loop': (259, 275), 'latch': (349, 354), '_offset':100}
    coords['1l5a'] = {'floor': (265, 271), 'activesite': (124, 131), 'latch_loop': (337, 361), 'B1': (21,30), 'Nlobe':(0,173), 'A5': (159, 173), 'A6': (173, 190), 'Clobe':(173,424), 'floor_loop': (260, 276), 'latch': (347, 356),  '_offset':"FIRST"}
    coords['4jn3'] = {'floor': (298, 304), 'activesite': (156, 163), 'latch_loop': (367, 392), 'B1': (32,40), 'Nlobe':(0,207), 'A5': (192, 207), 'A6': (207, 221), 'Clobe':(207,450), 'floor_loop': (294, 310), 'latch': (378, 386),  '_offset':"FIRST"}
    refcoords = coords[ref]

    allresults = {}
    if not ncpus:
        ncpus = mp.cpu_count()
    pool = mp.Pool(ncpus)
    for i,x in enumerate(filenames):
        if os.path.exists(x):
            allresults[x]=pool.apply_async(processclustgbks_ASv5, (x,), dict(refpdb=refpdb,refseq=refseq,mergetol=mergetol,alignmethod=alignmethod,refcoords=refcoords,savealignments=savealignments,alignfolder=alignfolder))
    pool.close()
    pool.join()
    #collect results
    results = [item for sublist in allresults.values() for item in sublist.get() if sublist.successful()]
    if db:
        #flatten results and write to db
        sqlmethods.writetodb(results,db)
    return results

global Cdomextend # number of basepairs to extend Cdomain model
Cdomextend=450
global log
log = setlog.init(toconsole=True)