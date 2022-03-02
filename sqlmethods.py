#!/usr/bin/env python
# 2019 Mohammad Alanjary
# Wageningen University and Research (WUR)
# Bioinformatics department, Lab of Prof. Marnix Medema
#
# License: A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import sqlite3 as sql
import setlog

def writetodb(results,db):
    #write results to sql database
    conn = sql.connect(db)
    csr = conn.cursor()
    allkeys = set()
    for rec in results:
        allkeys |= set(rec.keys())
    allkeys.remove("modid")
    #Make table if doesnt exist:
    try:
        sqlcmd2= 'CREATE TABLE Modules ("modid" text PRIMARY KEY,'+", ".join(['"%s" text'%x for x in list(allkeys)])+")"
        sqlcmd='CREATE TABLE Modules ("modid" text PRIMARY KEY, "upid" text, "type" text, "seq" text, "modnum" text, "extraDomains" text, "Edom" text, "downid" text,\
         "domstring" text, "domseqs" text, "domlocs" text, "domlist" text, "domaincount" text, "clusterid" text, "C1type" text, "C1Seq#Nlobe" text, "C1Seq#latch_loop" text,\
          "C1Seq#latch" text, "C1Seq#floor_loop" text, "C1Seq#floor" text, "C1Seq#end" text, "C1Seq#Clobe" text, "C1Seq#B1" text, "C1Seq#activesite" text, "C1Seq#A6" text,\
           "C1Seq#A5" text, "C1PSeq#Nlobe" text, "C1PSeq#latch_loop" text, "C1PSeq#latch" text, "C1PSeq#floor_loop" text, "C1PSeq#floor" text, "C1PSeq#end" text, "C1PSeq#Clobe" text,\
            "C1PSeq#B1" text, "C1PSeq#activesite" text, "C1PSeq#A6" text, "C1PSeq#A5" text, "C1loc#Nlobe" text, "C1loc#latch_loop" text, "C1loc#latch" text, "C1loc#floor_loop" text,\
             "C1loc#floor" text, "C1loc#end" text, "C1loc#Clobe" text, "C1loc#B1" text, "C1loc#activesite" text, "C1loc#A6" text, "C1loc#A5" text, "C2type" text, "C2Seq#Nlobe" text,\
              "C2Seq#latch_loop" text, "C2Seq#latch" text, "C2Seq#floor_loop" text, "C2Seq#floor" text, "C2Seq#end" text, "C2Seq#Clobe" text, "C2Seq#B1" text, "C2Seq#activesite" text,\
               "C2Seq#A6" text, "C2Seq#A5" text, "C2PSeq#Nlobe" text, "C2PSeq#latch_loop" text, "C2PSeq#latch" text, "C2PSeq#floor_loop" text, "C2PSeq#floor" text, "C2PSeq#end" text,\
                "C2PSeq#Clobe" text, "C2PSeq#B1" text, "C2PSeq#activesite" text, "C2PSeq#A6" text, "C2PSeq#A5" text, "C2loc#Nlobe" text, "C2loc#latch_loop" text, "C2loc#latch" text,\
                 "C2loc#floor_loop" text, "C2loc#floor" text, "C2loc#end" text, "C2loc#Clobe" text, "C2loc#B1" text, "C2loc#activesite" text, "C2loc#A6" text, "C2loc#A5" text, \
                 "C2loc#acceptNlobe1" text, "C1loc#acceptNlobe1" text, "C2PSeq#acceptNlobe1" text, "C1PSeq#acceptNlobe1" text, "C2Seq#acceptNlobe1" text, "C1Seq#acceptNlobe1" text,\
                 "C2loc#acceptNlobe2" text, "C1loc#acceptNlobe2" text, "C2PSeq#acceptNlobe2" text, "C1PSeq#acceptNlobe2" text, "C2Seq#acceptNlobe2" text, "C1Seq#acceptNlobe2" text,\
                 "C2loc#acceptClobe1" text, "C1loc#acceptClobe1" text, "C2PSeq#acceptClobe1" text, "C1PSeq#acceptClobe1" text, "C2Seq#acceptClobe1" text, "C1Seq#acceptClobe1" text,\
                 "C2loc#acceptClobe2" text, "C1loc#acceptClobe2" text, "C2PSeq#acceptClobe2" text, "C1PSeq#acceptClobe2" text, "C2Seq#acceptClobe2" text, "C1Seq#acceptClobe2" text,\
                  "AdomSpecMethod" text, "AdomSpec" text, "accession" text)'
        log.debug(sqlcmd2)
        csr.execute(sqlcmd2)
        log.info("Creating Modules Table...")
    except sql.OperationalError as ex:
        log.info("Error: %s - Modules table exists. Adding to database..."%ex)
    for result in results:
        cols = '"'+'", "'.join(result.keys())+'"'
        sqlcmd="INSERT INTO Modules (%s) VALUES (%s)"%(cols,', '.join('?'*len(result)))
        csr.execute(sqlcmd,[str(x) for x in result.values()])
    conn.commit()

global log
log = setlog.init(toconsole=True)