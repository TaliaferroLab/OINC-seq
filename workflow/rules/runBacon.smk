"""
Rules for running bacon.py
For usage, include this in your workflow.
"""

rule runBacon:
    """runs Bacon with desired sample comparisons"""
    input:
        bacon=config["bacon"]["bacon"],
        pigpen=expand("PIGPEN_alignments/{outDir}/{sample}.pigpen.txt", sample=config["samples"], outDir =config["pigpen"]["outputDir"]),
    params:
        SC=expand("{sampconds}", sampconds=config["bacon"]["sampconds"]),
        MR=" --minreads " + config["bacon"]["minreads"] if config["bacon"]["minreads"] else "",
        CA=" --conditionA " + config["bacon"]["conditionA"] if config["bacon"]["conditionA"] else "",
        CB=" --conditionB " + config["bacon"]["conditionB"] if config["bacon"]["conditionB"] else "",
        OT=" " + config["bacon"]["tags"] if config["bacon"]["tags"] else "",
        OD=expand("PIGPEN_alignments/{outDir}/{output}", output=config["bacon"]["output"], outDir =config["pigpen"]["outputDir"]),
    output:
        expand("PIGPEN_alignments/{outDir}/{output}.bacon.txt", output=config["bacon"]["output"], outDir =config["pigpen"]["outputDir"]),
    threads: config["threads"]
    run:
        for sampcond,output in zip(params.SC,params.OD):
            shell(
                "python -u {input.bacon} --sampconds {sampcond} "
                "--output {output}{params.MR}{params.CA}{params.CB}{params.OT}"
            )
            
