import os
# from jinja2 import Template
# import rpy2.robjects as robjects 
import pypeliner.commandline 


def sample_level_report(mutations_per_cell, summary, 
                      snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
                      snv_cell_counts, snv_alt_counts, rearranegementtype_distribution, chromosome_types,
                      BAFplot, CNplot, datatype_summary, html_file, out_dir, sample_id):

	# cwd = os.path.dirname(os.path.abspath(__file__))
	# print(CNplot, html_file)
	# rcode = "rmarkdown::render('{}', output_file = '{}', params=list(cn_plot='{}'))".format("/juno/work/shah/abramsd/CODE/single_cell_pipeline/single_cell/workflows/qc/scripts/report2.Rmd", html_file, CNplot)
	# cmd = ["R", "-e",  rcode]
	# pypeliner.commandline.execute(*cmd, docker_image="rocker/tidyverse")

	rcode = "rmarkdown::render('/juno/work/shah/abramsd/CODE/single_cell_pipeline/single_cell/workflows/qc/scripts/report2.Rmd',output_file = '{}',  params=list( sample_id='{}', mutations_per_cell_png='{}', summary_csv='{}', snvs_high_impact_csv='{}', snvs_all_csv='{}', trinuc_csv='{}', snv_adjacent_distance_png='{}', snv_genome_count_png='{}', snv_cell_counts_png='{}', snv_alt_counts_png='{}', rearranegementtype_distribution_png='{}', chromosome_types_png='{}', BAFplot_png='{}', cn_plot_png='{}', datatype_summary_csv='{}'))".format(
					html_file, sample_id, mutations_per_cell, summary, 
					snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
					snv_cell_counts, snv_alt_counts, rearranegementtype_distribution, chromosome_types,
					BAFplot, CNplot, datatype_summary)


	cmd = ["R", "-e",  rcode]
	pypeliner.commandline.execute(*cmd, docker_image="rocker/tidyverse")
	# rcode = '''Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")
	# rmarkdown::render("scripts/report.Rmd", output_file = {}, params=list(args = myarg))'''
	# rcode = rcode.format(html_file)
	
	# robjects.r(rcode)
    # pypeliner.commandline.execute(
    #     'fastqc',
    #     '--outdir=' + temp_dir,
    # cmd = ['samtools', 'index', bamfile, indexfile]


def sample_level_report2(report, html_file, mutations_per_cell, summary, 
                      snvs_high_impact, snvs_all, trinuc, snv_adjacent_distance, snv_genome_count, 
                      snv_cell_counts, snv_alt_counts, rearranegementtype_distribution, chromosome_types,
                      BAFplot, CNplot, datatype_summary):

	cwd = os.path.dirname(os.path.abspath(__file__))
	template = Template(open(os.path.join(cwd, "scripts/sample_level_report.html"),"r").read())
	report = os.path.join(html_file)

	html = template.render(rearranegementtype_distribution=rearranegementtype_distribution, 
						   chromosome_types=chromosome_types, 
						   mutations_per_cell = mutations_per_cell, 
						   summary=summary,  snvs_high_impact = snvs_high_impact, snvs_all =snvs_all, 
						   trinuc =trinuc, snv_adjacent_distance=snv_adjacent_distance, snv_genome_count=snv_genome_count, 
                      	   snv_cell_counts=snv_cell_counts, snv_alt_counts=snv_alt_counts, BAFplot=BAFplot, CNplot=CNplot, 
						   datatype_summary=datatype_summary)
	output = open(report,"w")
	output.write(html)
	output.close()


