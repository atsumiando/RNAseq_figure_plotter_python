#please type this code "conda install -c anaconda seaborn=0.9.0" to update seaborn to use rnaseq_figure_plotter software.



import os
import time
import argparse
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.mlab as mlab
from scipy import stats

version = "0.0.2"

#parameter settings
parser = argparse.ArgumentParser(usage='python rnaseq_figure_plotter.py -i input_file -t bar -o output_file -g gene_list_file  ... -c 5 -s 6', add_help=True,formatter_class=argparse.RawDescriptionHelpFormatter)
#required parameter
parser.add_argument('-i','--input',help = 'input file name', required=True)
parser.add_argument('-t', '--type',help = 'choose plot types (bar, box, density, dot, heatmap, histogram, line, scatter, or violin)',required=True)
#general optional parameter
parser.add_argument('-o', '--output', help = 'default output; output file name',default = "output")
parser.add_argument('-g', '--gene',help = 'file name of specific gene ID list; generate "output"_gene_selection.txt file')
parser.add_argument('-l', '--log',help = 'default None; calculate log value (log2; 2, log10; 10, loge; e)',default = '')
parser.add_argument('-lgn', '--log_number',help = 'default 0.000000001; add number to avoid -inf for log value',default = 0.000000001, type = float)
parser.add_argument('-x', '--xaxis',help = 'default samples; choose x-axis (gene, sample, or data)',default = "sample")
parser.add_argument('-y', '--yaxis',help = 'default data; choose y-axis (gene, sample, or data)',default = "data")
parser.add_argument('-z', '--zaxis',help = 'default gene; choose z-axis (gene, sample, or data)',default = "gene")
parser.add_argument('-c', '--color',help = 'default 1; choose color type (1-10)', default = 1, type = int)
parser.add_argument('-f', '--figure_save_format',help = 'default pdf; choose format of figures (eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, or tiff)', default = 'pdf')
#optional parameter for individual plot types
parser.add_argument('-s', '--style',help = 'default 1; choose style of figures (1-8)', default = 1, type = int)
parser.add_argument('-zs', '--zscore',help = 'default None; apply z-score transformation in heatmap. Z-score application in column or row  is --xaxis (column); 1, and --zaxis (row); 2)',default = [],type = int)
parser.add_argument('-cc', '--cluster_column',help = 'default None; apply column cluster function for heatmap (on; 1)',default = [])
parser.add_argument('-cr', '--cluster_row',help = 'default None; apply row cluster function for heatmap (on; 1)',default = [])
parser.add_argument('-sc', '--scatter_column',help = 'default None; type column of two samples for comparison in dot plot. Split samples by comma(,). (example "sample1,sample2")')
parser.add_argument('-sr', '--scatter_row',help = 'default None; type row of two genes for comparison in dot plot. Split genes by comma(,).(example "geneA,geneB")')
args = parser.parse_args()



#set parameter name
input_file_name = args.input
output_file_name = args.output
list_file_name = args.gene
x = args.xaxis
y = args.yaxis
z = args.zaxis
save_format = args.figure_save_format



#select gene expression data from gene list and generate gene_selection_file and data_frame_change_file
def selection_file():
        with open(input_file_name, 'r') as sourcefile:
                list = sourcefile.read().splitlines()
                inputs, record = [], []
                for index in range(len(list)):
                        record = list[index].split('\t')
                        inputs.append(record)
	with open(list_file_name, 'r') as sourcefile2:
                list2 = sourcefile2.read().splitlines()
                inputs2, record2 = [], []
                for index2 in range(len(list2)):
                        record2 = list2[index2].split('\t')
                        inputs2.append(record2)
        with open(output_file_name+'_gene_selection.txt', 'w' ) as sinkfile:
                for index in range(len(list)):
                	for index2 in range(len(list2)):
				if inputs[index][0] == inputs2[index2][0]:
                        		data = inputs[index]
                        		data_str= [str(x) for x in data]
                        		data_final ='\t'.join(str(e) for e in data_str)
                        		sinkfile.write(data_final +"\n")
	with open(output_file_name+'_gene_selection.txt', 'r' ) as sourcefile3:
                list3 = sourcefile3.read().splitlines()
                inputs3, record3 = [], []
                for index3 in range(len(list3)):
                        record3 = list3[index3].split('\t')
                        inputs3.append(record3)
	with open('data_frame_change_file', 'w' ) as sinkfile:
		for index3 in range(0,len(list3)):
			for i in range (1,len(inputs3[0])):
				if args.log =='2':
					data1 = inputs[0][i]
					data2 = np.log2(float(inputs3[index3][i])+args.log_number)
					data3 = inputs3[index3][0]
				elif args.log =='10':
					data1 = inputs[0][i]
					data2 = np.log10(float(inputs3[index3][i])+args.log_number)
					data3 = inputs3[index3][0]
				elif args.log =='e':
					data1 = inputs[0][i]
					data2 = np.log(float(inputs3[index3][i])+args.log_number)
					data3 = inputs3[index3][0]
				else:
					data1 = inputs[0][i]
					data2 = inputs3[index3][i]
					data3 = inputs3[index3][0]
				sinkfile.write(data1+"\t"+ str(data2)+"\t"+ data3 +"\n")



#generate data_frame_change_file
def no_selection_file():
	with open(input_file_name, 'r') as sourcefile:
                list = sourcefile.read().splitlines()
                inputs, record = [], []
                for index in range(len(list)):
                        record = list[index].split('\t')
                        inputs.append(record)
	with open('data_frame_change_file', 'w' ) as sinkfile:
		for index in range(1,len(list)):
			for i in range (1,len(inputs[0])):
				if args.log =='2':
					data1 = inputs[0][i]
					data2 = np.log2(float(inputs[index][i])+args.log_number)
					data3 = inputs[index][0]
				elif args.log =='10':
					data1 = inputs[0][i]
					data2 = np.log10(float(inputs[index][i])+args.log_number)
					data3 = inputs[index][0]
				elif args.log =='e':
					data1 = inputs[0][i]
					data2 = np.log(float(inputs[index][i])+args.log_number)
					data3 = inputs[index][0]
				else:
					data1 = inputs[0][i]
					data2 = inputs[index][i]
					data3 = inputs[index][0]
				sinkfile.write(data1+"\t"+ str(data2)+"\t"+ data3+"\n")



#select gene expression data from gene list and generate gene_selection_file, remove same value among all samples to generate data_frame_change_file
def selection_file_rm_duplicate():
    with open(input_file_name, 'r') as sourcefile:
        list = sourcefile.read().splitlines()
        inputs, record = [], []
        for index in range(len(list)):
            record = list[index].split('\t')
            inputs.append(record)
    with open(list_file_name, 'r') as sourcefile2:
        list2 = sourcefile2.read().splitlines()
        inputs2, record2 = [], []
        for index2 in range(len(list2)):
            record2 = list2[index2].split('\t')
            inputs2.append(record2)
        with open(output_file_name+'_gene_selection.txt', 'w' ) as sinkfile:
            for index in range(len(list)):
                for index2 in range(len(list2)):
                    if inputs[index][0] == inputs2[index2][0]:
                        data = inputs[index]
                        data_str= [str(x) for x in data]
                        data_final ='\t'.join(str(e) for e in data_str)
                        sinkfile.write(data_final +"\n")
    with open(output_file_name+'_gene_selection.txt', 'r' ) as sourcefile3:
        list3 = sourcefile3.read().splitlines()
        inputs3, record3 = [], []
        for index3 in range(len(list3)):
            record3 = list3[index3].split('\t')
            inputs3.append(record3)
    with open('duplicate_remove_file', 'w' ) as sinkfile2:
        for index3 in range(len(list3)):
            arrays= inputs3[index3]
            if not all(arrays[i] == arrays[i+1] for i in range(1,len(arrays)-1)):
                data_str= [str(x) for x in arrays]
                data_final ='\t'.join(str(e) for e in data_str)
                sinkfile2.write(data_final +"\n")
    with open('duplicate_remove_file', 'r') as sourcefile4:
        list4 = sourcefile4.read().splitlines()
        inputs4, record4 = [], []
        for index4 in range(len(list4)):
            record4 = list4[index4].split('\t')
            inputs4.append(record4)
    with open('data_frame_change_file', 'w' ) as sinkfile3:
        for index4 in range(0,len(list4)):
            for i in range (1,len(inputs4[0])):
                if args.log =='2':
                    data1 = inputs[0][i]
                    data2 = np.log2(float(inputs4[index4][i])+args.log_number)
                    data3 = inputs4[index4][0]
                elif args.log =='10':
                    data1 = inputs[0][i]
                    data2 = np.log10(float(inputs4[index4][i])+args.log_number)
                    data3 = inputs4[index4][0]
                elif args.log =='e':
                    data1 = inputs[0][i]
                    data2 = np.log(float(inputs4[index4][i])+args.log_number)
                    data3 = inputs4[index4][0]
                else:
                    data1 = inputs[0][i]
                    data2 = inputs4[index4][i]
                    data3 = inputs4[index4][0]
                sinkfile3.write(data1+"\t"+ str(data2)+"\t"+ data3 +"\n")
    os.remove('duplicate_remove_file') #remove dataframe file


#remove same value among all samples to generate data_frame_change_file
def no_selection_file_rm_duplicate():
    with open(input_file_name, 'r') as sourcefile:
        list = sourcefile.read().splitlines()
        inputs, record = [], []
        for index in range(len(list)):
            record = list[index].split('\t')
            inputs.append(record)
    with open('duplicate_remove_file', 'w' ) as sinkfile:
        for index in range(len(list)):
            arrays= inputs[index]
            if not all(arrays[i] == arrays[i+1] for i in range(1,len(arrays)-1)):
                data_str= [str(x) for x in arrays]
                data_final ='\t'.join(str(e) for e in data_str)
                sinkfile.write(data_final +"\n")
    with open('duplicate_remove_file', 'r') as sourcefile2:
        list2 = sourcefile2.read().splitlines()
        inputs2, record2 = [], []
        for index2 in range(len(list2)):
            record2 = list2[index2].split('\t')
            inputs2.append(record2)
    with open('data_frame_change_file', 'w' ) as sinkfile2:
        for index2 in range(1,len(list2)):
            for i in range (1,len(inputs2[0])):
                if args.log =='2':
                    data1 = inputs2[0][i]
                    data2 = np.log2(float(inputs2[index2][i])+args.log_number)
                    data3 = inputs2[index2][0]
                elif args.log =='10':
                    data1 = inputs2[0][i]
                    data2 = np.log10(float(inputs2[index2][i])+args.log_number)
                    data3 = inputs2[index2][0]
                elif args.log =='e':
                    data1 = inputs2[0][i]
                    data2 = np.log(float(inputs2[index2][i])+args.log_number)
                    data3 = inputs2[index2][0]
                else:
                    data1 = inputs[0][i]
                    data2 = inputs2[index2][i]
                    data3 = inputs2[index2][0]
                sinkfile2.write(data1+"\t"+ str(data2)+"\t"+ data3+"\n")
    os.remove('duplicate_remove_file') #remove dataframe file



#select types of data_frame_change_file
def selection_or_not():
    if bool(list_file_name):
        if bool(args.cluster_column):
            return 	selection_file_rm_duplicate()
        if bool(args.cluster_row):
            return  selection_file_rm_duplicate()
        else:
            return selection_file()
    else:
        if bool(args.cluster_column):
            return no_selection_file_rm_duplicate()
        if bool(args.cluster_row):
            return no_selection_file_rm_duplicate()
        else:
            return no_selection_file()



#select types of graph (--type)
def graph_types():
    if args.type == "scatter":
        if bool(args.scatter_row):
            return scatter_row()
        elif bool(args.scatter_column):
            return scatter_column()
        else:
            print("error; -sc or -sr are not provided")
    elif args.type == "heatmap":
        if bool(args.zscore):
            return heatmap_plot_zscore()
        else:
            return heatmap_plot()
    elif args.type == "line":
        return line_plot()
    elif args.type == "box":
        return box_plot()
    elif args.type == "bar":
        return bar_plot()
    elif args.type == "dot":
        return dot_plot()
    elif args.type == "violin":
        return violin_plot()
    elif args.type == "density":
        return density_plot()
    elif args.type == "histogram":
        return histogram_plot()
    else:
        print("error; graph type(-t, --type) are not provided")



#select types of color palette (--color)
def color_change():
	if args.color == 1:
		letter = "RdBu_r"
	elif args.color == 2:
		letter = "Reds"
	elif args.color == 3:
		letter = "Blues"
	elif args.color == 4:
		letter = "RdYlBu_r"
	elif args.color == 5:
		letter = "RdGy_r"
	elif args.color == 6:
		letter = "Paired"
	elif args.color == 7:
		letter = "cubehelix"
	elif args.color == 8:
		letter = "muted"
	elif args.color == 9:
		letter = "hls"
	else:
		letter = "Set2"
	return letter  



#select types of style; background (--style)
def style_change():
	if args.style <= 2:
		return ("whitegrid")
	elif args.style  <= 4 :
		return ("white")
	elif args.style <=6:
		return ("darkgrid") 
	else:
		return ("dark")



#select types of style; size (--style)
def style_change2():
	mod = args.style%2
	if mod == 1 :
		letter = "paper"
	else:
		letter = "talk"
	return letter



#bar plot
def bar_plot():
	figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
	figure_df.columns = ['sample', 'data', 'gene']
	sns.set_style(style_change())
	sns.set_context(style_change2())
	figure = sns.barplot(x=x, y=y, hue=z, data=figure_df, palette= color_change())
	figure.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)
	plt.savefig(output_file_name +"_barplot." +save_format, bbox_inches='tight',dpi=100)
	plt.close()



#box plot
def box_plot():
	figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
	figure_df.columns = ['sample', 'data', 'gene']
	sns.set_style(style_change())
	sns.set_context(style_change2())
	figure = sns.boxplot(x=x, y=y, data=figure_df, palette= color_change())
	plt.savefig(output_file_name +"_boxplot." +save_format, bbox_inches='tight',dpi=100)
	plt.close()


#density plot
def density_plot():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    samples = list(figure_df[x].unique())
    sns.set_style(style_change())
    sns.set_context(style_change2())
    for sample in samples:
        subset = figure_df[figure_df[x] == sample]
        with sns.color_palette(color_change(),len(samples)):
            figure = sns.distplot(subset['data'], hist = False, kde = True, norm_hist = True, label = sample )
    figure.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)
    plt.ylabel('density')
    plt.savefig(output_file_name +"_densityplot." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#dot plot
def dot_plot():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    sns.set_style(style_change())
    sns.set_context(style_change2())
    figure = sns.scatterplot(x=x, y=y, hue=z, data=figure_df, palette= color_change())
    figure.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)
    plt.savefig(output_file_name +"_dotplot." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#heatmap without z_score
def heatmap_plot():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    figure_df= figure_df.pivot(z, x, y)
    figure = sns.clustermap(figure_df, cmap= color_change(), col_cluster= bool(args.cluster_column), row_cluster=bool(args.cluster_row))
    plt.savefig(output_file_name +"_heatmap." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#heatmap with z_score
def heatmap_plot_zscore():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    figure_df= figure_df.pivot(z, x, y)
    figure = sns.clustermap(figure_df, cmap= color_change(), z_score = args.zscore, col_cluster= bool(args.cluster_column), row_cluster=bool(args.cluster_row))
    plt.savefig(output_file_name +"_heatmap." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#histogram plot
def histogram_plot():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    samples = list(figure_df[x].unique())
    sns.set_style(style_change())
    sns.set_context(style_change2())
    for sample in samples:
        subset = figure_df[figure_df[x] == sample]
        with sns.color_palette(color_change(),len(samples)):
            figure = sns.distplot(subset['data'], hist = True, kde = False, label = sample )
    figure.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)
    plt.ylabel('count')
    plt.savefig(output_file_name +"_histogram." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#line plot
def line_plot():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    sns.set_style(style_change())
    sns.set_context(style_change2())
    figure = sns.lineplot(x=x, y=y, hue=z, data=figure_df, palette= color_change())
    figure.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)
    plt.savefig(output_file_name +"_lineplot." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#scatter plot figure of column
def scatter_column():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    figure_df= figure_df.pivot(z, x, y)
    pair_c = []
    pair_c = args.scatter_column.split(",")
    col1 = figure_df.loc[:,pair_c[0]]
    col2 = figure_df.loc[:,pair_c[1]]
    sns.set_style(style_change())
    sns.set_context(style_change2())
    figure = sns.scatterplot(x=col1, y=col2)
    plt.savefig(output_file_name+"_" +pair_c[0]+"_"+pair_c[1] +"_scatterplot." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#scatter plot figure of row
def scatter_row():
    figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
    figure_df.columns = ['sample', 'data', 'gene']
    figure_df= figure_df.pivot(z, x, y)
    pair_r = []
    pair_r = args.scatter_row.split(",")
    row1 = figure_df.loc[pair_r[0]]
    row2 = figure_df.loc[pair_r[1]]
    sns.set_style(style_change())
    sns.set_context(style_change2())
    figure = sns.scatterplot(x=row1, y=row2)
    plt.savefig(output_file_name+"_"+pair_r[0]+"_"+pair_r[1]  +"_scatterplot." +save_format, bbox_inches='tight',dpi=100)
    plt.close()



#violin plot
def violin_plot():
	figure_df = pd.read_csv("data_frame_change_file", sep='\t', header=None)
	figure_df.columns = ['sample', 'data', 'gene']
	sns.set_style(style_change())
	sns.set_context(style_change2())
	figure = sns.violinplot(x=x, y=y, data=figure_df, palette= color_change())
	plt.savefig(output_file_name +"_violinplot." +save_format, bbox_inches='tight',dpi=100)
	plt.close()



#start operation
start_time = time.time()
print("You are using rnaseq_figure_plotter version; %s" % version)
selection_or_not()
print("Complete dataframe generation! Prepareing for plotting!")
print("Dataframe generation time; %s seconds" % round((time.time() - start_time),2))
graph_types()
os.remove('data_frame_change_file') #remove dataframe file
print("Plotting now!")
print("Total operation time; %s seconds" % round((time.time() - start_time),2))
