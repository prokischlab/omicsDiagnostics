import pathlib
import wbuild
config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)

configfile: "wbuild.yaml"
include: config['wBuildPath'] + "/wBuild.snakefile"
htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

rule all:
    input: rules.Index.output, "Output/html/readme.html"
    output: touch("Output/all.done")

rule outrider:
    input: config["PROC_DATA"] + "/outrider/OUTRIDER_results.rds"
    output: touch("Output/outrider.done")
        
rule limma:
    input: config["PROC_DATA"] + "/limma/LIMMA_results.rds"
    output: touch("Output/limma.done")
    
rule protrider:
    input: config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"
    output: touch("Output/protrider.done")