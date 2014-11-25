#!/usr/bin/perl


###########################################################################
## ## A script to generate an HTML main document with CAP-miRSeq results
##
## 
##
## Parameters:
## <run_info.txt> - The CAP-miRSeq run_info.txt config file
## <output file> - ouput HTML file name
## <sample_summary.txt> - table of per-sample statistics
## <SNVs called (0/1)> - Boolean to indicate whether SNVs were called in this analysis
## <trim adapter(0/1)> - Boolean to indicate whether adapters were trimmed in this analysis
## <diff expression analysis> - Optional colon seperated list of differential expression labels
##
############################################################################


use strict;
use warnings;

open RUN_INFO, "$ARGV[0]" or die "opening file: $ARGV[0]";
open OUT, ">", "$ARGV[1]" or die "opening file: $ARGV[1]";
open STATS, "$ARGV[2]" or die "opening file: $ARGV[2]";
#open MIRNA, "$ARGV[3]" or die "opening file: $ARGV[3]"; 

# read config variables
#while(<RUN_INFO>){
#	my $row = $_;
#	chomp $row;
#	my @line = split("=",$row);
#	$config{$line[0]} = $line[1];
#}

#my @samples = split(";",$config{"SAMPLENAMES"});



print OUT <<EOM;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>Onyx | miRNA Differential Expression Analysis Pipeline</title>
<style type="text/css">
html {
	height: 100%;
}
body {
	width: 100%;
	height: 100%;
	background: #003333;
	background-repeat: repeat;
	text-align: left;
	margin: 0px 0px 0px 0px;
	padding: 0px 0px 0px 0px;
	font-family:verdana, Arial, Helvetica, sans-serif;
	font-weight: normal;
	color: #333;
	font-size: 12px;
}
h1,h2,h3 {
	font-family:Georgia, "Times New Roman", Times, serif;
}
ul {
	margin-bottom: 10px;
}
ul li {
	margin-bottom: 5px;
}
table {
	border: none;
}
td {
	padding: 5px 15px 5px 15px;
	margin: 0px 0px 0px 0px;
}
tr {

}
tr.shade {
	background: #ddd;
}
.strong {
	font-weight: bold;
	background: #ddd;
}
table.summary {
	font-size: 12px;
}
#wrapper {
	margin: 0px auto 0px auto;
	padding: 0px 50px 10px 100px;
}
#header {
	/*width: 800px;*/
	height: 100px;
	margin: 0px auto 0px auto;
	padding: 0px 50px 0px 100px;
	/*background: #444;*/
	/*background: #336699;*/
	background: #506987;
	text-align: left;
	font-family:Georgia, "Times New Roman", Times, serif;
}
#header h1 {
	color: #fefefe;
	font-size: 50px;
	margin: 0px 0px 0px 0px;
	padding: 5px 5px 5px 0px;
}
#header h3 {
	color: #fefefe;
	font-size: 20px;
	margin: 0px 0px 0px 0px;
	padding: 0px 5px 5px 0px;
}
#footer {
	margin: 0px auto 0px auto;
	padding: 10px 50px 20px 100px;
	background: #444;
	color: #fefefe;
}
#footer a{
	color: #fefefe;
}
</style>
</head>
<body>
<div id="header">
<h1>Onyx <span style="font-size: 25px;">v1.1</span></h1>
<h3>miRNA Differential Expression Analysis Pipeline</h3>
</div>
<div id="wrapper">
<h2>Project Summary</h2>
<table>

EOM

my @date = localtime(time);
my $day = $date[3];
my $month = $date[4]+1;
my $year = $date[5]+1900;

my $finger_bin = `which finger` ; 

my @name ; 

if($finger_bin ne ""){
	chomp(my @finger = `finger $ENV{"USER"}`);
	my @name_row = split("Name:",$finger[0]);
	@name = split(";",$name_row[1]);
} else {
	my $whoami = `whoami`; 
	chomp($whoami);
	push(@name, $whoami);
} 


#print OUT '<tr class="shade"><td>Date</td><td>'."$month/$day/$year".'</td></tr>
#<tr><td>Number of Samples</td><td>'.scalar(@samples).'</td></tr>
#<tr class="shade"><td>Genome Build</td><td>'.$config{"GENOME_BUILD"}.'</td></tr>
#<tr><td>miRBase Version</td><td>'.$config{"MIRBASE_VERSION"}.'</td></tr>
#<tr class="shade"><td>Analysis Performed By</td></tr>';

print OUT <<ENDHTML;
</table>
<h2>Analysis Workflow</h2>
<img src="onyx_pipeline.csv" alt="Onyx Pipeline Summary" style="width:768px; height:576px;" />

<h2>Quality Control Reports</h2>
<ul>


<li><a href="doc/pretrim"><strong>FASTQC Reports (Before Trimming)</strong></a></li>

<li><a href="doc/posttrim"><strong>FASTQC Reports (After Trimming)</strong></a></li>
<li><a href="qc/other_rna"><strong>Quantification of Other RNA</strong></a></li>
</ul>

<h2>Normalization Quality</h2>
<img src="result/expression/QUA_NOR_plot.png" alt="Normalization Quality" style="width:800px; " />

<h2>Sample Summary</h2>
<table class="summary">

<h2>MDS Plot</h2>
<img src="result/expression/MDS_Plot.png" alt="MDS Plot " style="width:800px; " />

ENDHTML

# print sample summary table
my $header = 1;
my $shade = 1;
while(<STATS>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my $class = "";
	if($header){
		$class = "class=\"strong\"";
		$header = 0;
	}
	if($shade){
		print OUT "<tr $class class=\"shade\">";
		$shade = 0;
	}else{
		print OUT "<tr>";
		$shade = 1;
	}
	foreach my $value (@line){
		print OUT "<td>".$value."</td>";
	}
	print OUT "</tr>\n";

}

print OUT <<ENDHTML;
</table>
<h2>miRNas with PValue > 0.05 , logFC > |1| and FDR < 0.01 </h2>
<img src="barplot_plot.png" alt="miRNAs bar " style="width:800px; " />

<h2>Heatmap of miRNas with PValue > 0.05, logFC > |1| and FDR < 0.01</h2>
<img src="heatmap_plot.png" alt="miRNAs bar " style="width:800px; " />

ENDHTML

# my $header2 = 1;
# my $shade2 = 1;
# while(<miRNAs_dif_expressed>){
# 	my $row2 = $_;
# 	chomp $row2;
# 	my @line2 = split("\t",$row2);
# 	my $class2 = "";
# 	if($header2){
# 		$class2 = "class=\"strong\"";
# 		$header2 = 0;
# 	}
# 	if($shade2){
# 		print OUT "<tr $class2 class=\"shade\">";
# 		$shade2 = 0;
# 	}else{
# 		print OUT "<tr>";
# 		$shade2 = 1;
# 	}
# 	foreach my $value2 (@line2){
# 		print OUT "<td>".$value2."</td>";
# 	}
# 	print OUT "</tr>\n";

# }


print OUT <<ENDHTML;
</table>
<h2>Result Reports</h2>
<ul>
<li><a href="result/expression"><strong>Expression Reports (Merged With All Samples)</strong></a></li>
<ul>
<li><a href="result/expression/differential_expression.xls">Raw expression.xls</a> - Raw Expression Results</li>
<li><a href="result/expression/differential_expression.xls">Differential_expression.xls</a> - Differential Expression Results</li>
<li><a href="result/expression/dif_exp_sum.xls">Differential expression of miRNAs significant.xls</a> - Differential Expression Results of miRNas with PValue > 0.05 , logFC > |1| and FDR < 0.01</li>
</ul>


<h2>FDR Volcano Plot</h2>
<img src="result/expression/FDR.volcano_plot.png" alt="FDR Volcano Plot" style="width:800px; " />

<h2>P-Value Vocano Plot</h2>
<img src="result/expression/PValue_plot.png" alt="PValue Volcano Plot" style="width:800px; " />

<h2>P-Value distribution</h2>
<img src="result/expression/PValue_dis_plot.png" alt="PValue distribution" style="width:800px; " />

<h2>SMEAR Plot</h2>
<img src="result/expression/SMEAR_plot.png" alt="SMEAR Plot" style="width:800px; " />

ENDHTML

print OUT <<ENDHTML;
<h2>View Alignments in IGV</h2>
<ul>
<li><a href="igv/igv_session.xml"><strong>igv_session.xml</strong></a> - Open this session file in <a href="http://www.broadinstitute.org/software/igv/download">IGV</a> to view aligned reads.</li>
</ul>
<h2>What Next?</h2>
<ul>
<li>Differentially expressed miRNAs among conditions</li>
<li>Integration with gene expression data</li>
<li>miRNA target identification</li>
<li>Pathway analysis (miRNA target canonical pathway or network analysis)</li>
<li>Variant/mutation identification in miRNA coding region</li>
</ul>
<h2>Useful Links</h2>
<ul>
<li><a href="http://bioinformaticstools.mayo.edu" target="_blank">CAP-miRseq</a></li>
<li><a href="http://www.mirbase.org" target="_blank">miRBase</a></li>
<li><a href="https://www.mdc-berlin.de/8551903/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep" target="_blank">miRDeep2</a></li>
<li><a href="http://www.broadinstitute.org/software/igv/download" target="_blank">Integrative Genomics Viewer (IGV)</a></li>
<li><a href="http://www.bioconductor.org/packages/2.13/bioc/html/edgeR.html" target="_blank">edgeR</a></li>

</ul>

</div>
<div id="footer">
<h2>Onyx Pipeline 2014.</h2>
</div>
</body>
</html>

ENDHTML


close RUN_INFO;
close OUT;
close STATS;
#close MIRNA;

