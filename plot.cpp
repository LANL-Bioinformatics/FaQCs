#include "FaQCs.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <unistd.h> // Needed for unlink

using namespace std;

void write_matrix(const string &m_filename, const matrix<size_t> &m_data);
void write_quality_histogram(const string &m_filename, const vector<size_t> &m_read_hist,
	const vector<size_t> &m_base_hist);
void write_base_content(const string &m_filename, const vector<NucleotideCount> &m_data);
void write_length_histogram(const string &m_filename, const vector<size_t> &m_hist);
void write_kmer_histogram(const string &m_filename, const MAP<size_t, size_t> &m_hist);
void write_rarefaction(const string &m_filename, const vector<Rarefaction> &m_rarefaction);

void plot(const PlotInfo &m_info, vector<size_t> &m_filter_stats, const Options &m_opt)
{
	// Use 'R' to produce a PDF document with graphs and information about
	// the reads we've just processed. We will need to feed the data to 'R'
	// using temporary files.
	//
	// TODO: Remove the need for temporary files by embedding the data directly in the
	// R script file. This will avoid the need for disk access (a small potential speed up)
	// and decrease the attack surface (better computer security).
	
	///////////////////////////////////////////////////////////////////////
	const string pre_quality_matrix_file = m_opt.output_dir + PATH_SEPARATOR + "qa." +
		m_opt.prefix + ".quality.matrix";
	const string post_quality_matrix_file = m_opt.output_dir + PATH_SEPARATOR +
		m_opt.prefix + ".quality.matrix";
		
	write_matrix(pre_quality_matrix_file, m_info.pre_quality_matrix);
	write_matrix(post_quality_matrix_file, m_info.post_quality_matrix);
	
	///////////////////////////////////////////////////////////////////////
	
	const string pre_base_matrix_file = m_opt.output_dir + PATH_SEPARATOR + "qa." +
		m_opt.prefix + ".base.matrix";
	const string post_base_matrix_file = m_opt.output_dir + PATH_SEPARATOR +
		m_opt.prefix + ".base.matrix";
	
	write_matrix(pre_base_matrix_file, m_info.pre_base_matrix);
	write_matrix(post_base_matrix_file, m_info.post_base_matrix);
	
	///////////////////////////////////////////////////////////////////////
	
	const string pre_quality_histogram_file = m_opt.output_dir + PATH_SEPARATOR +
		"qa." + m_opt.prefix + ".for_qual_histogram.txt";
	const string post_quality_histogram_file = m_opt.output_dir + PATH_SEPARATOR +
		m_opt.prefix + ".for_qual_histogram.txt";
		
	write_quality_histogram(pre_quality_histogram_file, 
		m_info.pre_read_quality_histogram, m_info.pre_base_quality_histogram);
	
	write_quality_histogram(post_quality_histogram_file,
		m_info.post_read_quality_histogram, m_info.post_base_quality_histogram);
	
	///////////////////////////////////////////////////////////////////////
	const string pre_nuc_composition_file = m_opt.output_dir + PATH_SEPARATOR +
		"qa." + m_opt.prefix + ".base_content.txt";
	const string post_nuc_composition_file = m_opt.output_dir + PATH_SEPARATOR +
		m_opt.prefix + ".base_content.txt";
		
	write_base_content(pre_nuc_composition_file, m_info.pre_nuc_composition);
	write_base_content(post_nuc_composition_file, m_info.post_nuc_composition);
	
	///////////////////////////////////////////////////////////////////////
	const string pre_length_histogram_file = m_opt.output_dir + PATH_SEPARATOR +
		"qa." + m_opt.prefix + ".length_count.txt";
	const string post_length_histogram_file = m_opt.output_dir + PATH_SEPARATOR +
		m_opt.prefix + ".length_count.txt";
		
	write_length_histogram(pre_length_histogram_file, m_info.pre_length_histogram);
	write_length_histogram(post_length_histogram_file, m_info.post_length_histogram);

	///////////////////////////////////////////////////////////////////////
	const string kmer_rarefaction_file = m_opt.output_dir + PATH_SEPARATOR + m_opt.prefix + ".Kmercount.txt";
	//const string kmer_files = m_opt.output_dir + PATH_SEPARATOR + m_opt.prefix + ".KmerFiles.txt";
	const string kmer_histogram_file = m_opt.output_dir + PATH_SEPARATOR + m_opt.prefix + ".kmerH.txt";

	// Don't test for m_info.kmer_rarefaction, since this flag will be false once we have
	// collected enough data to complete the requested rarefaction curve.
	if(m_info.kmer_frequency_histogram.empty() == false){

		write_kmer_histogram(kmer_histogram_file, m_info.kmer_frequency_histogram);
		write_rarefaction(kmer_rarefaction_file, m_info.kmer_rarefaction);
	}

	///////////////////////////////////////////////////////////////////////
	stringstream script;
	
	// Turn off R warnings (comment this line to debug R script problems)
	script << "options(warn=-1)\n";
	
	script << "if(file.exists(\"" << pre_length_histogram_file << "\")){\n";
	script << "  pdf(file = \"" << m_opt.plots_file << "\",width=15,height=7)\n";
	script << "}else{\n";
	script << "  pdf(file = \"" << m_opt.plots_file << "\",width=10,height=8)\n";
	script << "}\n";
	script << "def.par <- par(no.readonly = TRUE) # get default parameters\n";
	script << "\n";
	script << "#Summary\n";
	script << "par(family=\"mono\")\n";
	script << "SummaryStats<-readLines(\"" << m_opt.stats_file << "\")\n";
	script << "plot(0:1,0:1,type=\'n\',xlab=\"\",ylab=\"\",xaxt=\'n\',yaxt=\'n\',bty=\'n\')\n";
	script << "if (" << (m_opt.qc_only ? 1 : 0) << "){\n";
	script << "  for (i in 1:length(SummaryStats)){\n";
	script << "     text(0.05,1-0.04*(i-1),SummaryStats[i],adj=0,font=2,cex=1)\n";
	script << "  }\n";
	script << "}else{\n";
	script << "  if (" << m_filter_stats[FilterStat::PAIRED_READ_NUMBER] << " > 0) {\n";
	script << "      adjust <-14\n";
	script << "      abline(h=0.73,lty=2)\n";
	script << "  }else{\n";
	script << "      adjust<-11\n";
	script << "      abline(h=0.85,lty=2)\n";
	script << "  }\n";
	script << "  for (i in 1:length(SummaryStats)){\n";
	script << "     if (i>5 && i<adjust){\n";
	script << "       text(0.45,1-0.035*(i-6),SummaryStats[i],adj=0,font=2,cex=0.9)\n";
	script << "     }else if(i >=adjust){\n";
	script << "       text(0.05,1-0.035*(i-6),SummaryStats[i],adj=0,font=2,cex=0.9)\n";
	script << "     }else{\n";
	script << "       text(0.05,1-0.035*(i-1),SummaryStats[i],adj=0,font=2,cex=0.9)\n";
	script << "     }\n";
	script << "  }\n";
	script << "}\n";
	script << "#title(paste(\"" << m_opt.prefix << "\",\"QC report\"),sub = 'DOE Joint Genome Institute/Los Alamos National Laboratory', adj = 0.5, col.sub='darkblue',font.sub=2,cex.sub=0.8)\n";
	script << "title(\"QC stats\")\n";
	script << "par(def.par)#- reset to default\n";
	script << "\n";
	script << "#length histogram\n";
	script << "length_histogram <- function(length_count_file, xlab,ylab){\n";
	script << "  lengthfile<-read.table(file=length_count_file)\n";
	script << "  lengthList<-as.numeric(lengthfile$V1)\n";
	script << "  lengthCount<-as.numeric(lengthfile$V2)\n";
	script << "  lenAvg<-sum(lengthList * lengthCount)/sum(lengthCount)\n";
	script << "  lenStd<-sqrt(sum(((lengthList - lenAvg)**2)*lengthCount)/sum(lengthCount))\n";
	script << "  lenMax<-max(lengthList[lengthCount>0])\n";
	script << "  lenMin<-min(lengthList[lengthCount>0])\n";
	script << "  totalReads<-sum(lengthCount)\n";
	script << "  barplot(lengthCount/1000000,names.arg=lengthList,xlab=xlab,ylab=ylab,cex.names=0.8)\n";
	script << "  legend.txt<-c(paste(\"Mean\",sprintf (\"%.2f\",lenAvg),\"±\",sprintf (\"%.2f\",lenStd)),paste(\"Max\",lenMax),paste(\"Min\",lenMin))\n";
	script << "  legend('topleft',legend.txt,bty='n')\n";
	script << "  if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "      mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)\n";
	script << "  }\n";
	script << "  return(totalReads)\n";
	script << "}\n";
	script << "if(file.exists(\"" << pre_length_histogram_file << "\")){\n";
	script << "    par(mfrow=c(1,2),mar=c(5,6,4,2))\n";
	script << "    qa.readsCount<-length_histogram(\"" << pre_length_histogram_file << "\",\"Input Length\",\"Count (millions)\")\n";
	script << "    readsCount<-length_histogram(\"" << post_length_histogram_file << "\",\"Trimmed Length\",\"\")\n";
	script << "}else{\n";
	script << "    readsCount<-length_histogram(\"" << post_length_histogram_file << "\",\"Length\",\"Count (millions)\")\n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title(\"Reads Length Histogram\")\n";
	script << "\n";
	script << "#readGC plot\n";
	script << "\n";
	script << "readGC_plot <- function(base_content_file, totalReads, fig_x_start,fig_x_end,xlab,ylab,new){\n";
	script << "	baseP<-read.table(file=base_content_file)\n";
	script << "	Apercent<-baseP$V2[which(baseP$V1==\"A\")]\n";
	script << "	ApercentCount<-baseP$V3[which(baseP$V1==\"A\")]\n";
	script << "	Tpercent<-baseP$V2[which(baseP$V1==\"T\")]\n";
	script << "	TpercentCount<-baseP$V3[which(baseP$V1==\"T\")]\n";
	script << "	Cpercent<-baseP$V2[which(baseP$V1==\"C\")]\n";
	script << "	CpercentCount<-baseP$V3[which(baseP$V1==\"C\")]\n";
	script << "	Gpercent<-baseP$V2[which(baseP$V1==\"G\")]\n";
	script << "	GpercentCount<-baseP$V3[which(baseP$V1==\"G\")]\n";
	script << "	#Npercent<-baseP$V2[which(baseP$V1==\"N\")]\n";
	script << "	#NpercentCount<-baseP$V3[which(baseP$V1==\"N\")]\n";
	script << "	GCpercent<-baseP$V2[which(baseP$V1==\"GC\")]\n";
	script << "	GCpercentCount<-baseP$V3[which(baseP$V1==\"GC\")]\n";
	script << "	aAvg<-sum(Apercent * ApercentCount)/sum(ApercentCount)\n";
	script << "	aStd<-sqrt(sum(((Apercent - aAvg)**2)*ApercentCount)/sum(ApercentCount))\n";
	script << "	tAvg<-sum(Tpercent * TpercentCount)/sum(TpercentCount)\n";
	script << "	tStd<-sqrt(sum(((Tpercent - tAvg)**2)*TpercentCount)/sum(TpercentCount))\n";
	script << "	cAvg<-sum(Cpercent * CpercentCount)/sum(CpercentCount)\n";
	script << "	cStd<-sqrt(sum(((Cpercent - cAvg)**2)*CpercentCount)/sum(CpercentCount))\n";
	script << "	gAvg<-sum(Gpercent * GpercentCount)/sum(GpercentCount)\n";
	script << "	gStd<-sqrt(sum(((Gpercent - gAvg)**2)*GpercentCount)/sum(GpercentCount))\n";
	script << "	#nAvg<-sum(Npercent * NpercentCount)/sum(NpercentCount)\n";
	script << "	#nStd<-sqrt(sum(((Npercent - nAvg)**2)*NpercentCount)/sum(NpercentCount))\n";
	script << "	gcAvg<-sum(GCpercent * GCpercentCount)/sum(GCpercentCount)\n";
	script << "	gcStd<-sqrt(sum(((GCpercent - gcAvg)**2)*GCpercentCount)/sum(GCpercentCount))\n";
	script << "	GCaggregate<-tapply(GCpercentCount,list(cut(GCpercent,breaks=c(seq(0,100,1)))),FUN=sum)\n";
	script << "	Aaggregate<-tapply(ApercentCount,list(cut(Apercent,breaks=c(seq(0,100,1)))),FUN=sum)\n";
	script << "	Taggregate<-tapply(TpercentCount,list(cut(Tpercent,breaks=c(seq(0,100,1)))),FUN=sum)\n";
	script << "	Caggregate<-tapply(CpercentCount,list(cut(Cpercent,breaks=c(seq(0,100,1)))),FUN=sum)\n";
	script << "	Gaggregate<-tapply(GpercentCount,list(cut(Gpercent,breaks=c(seq(0,100,1)))),FUN=sum)\n";
	script << "\n";
	script << "	par(fig=c(fig_x_start,(fig_x_end-fig_x_start)*0.75+fig_x_start,0,1),mar=c(5,6,4,2),xpd=FALSE,cex.main=1.2,new=new)\n";
	script << "	plot(GCaggregate/1000000,xlim=c(0,100),type=\"h\",lwd=4,xlab=paste(xlab,\"GC (%)\"),ylab=ylab,lend=2)\n";
	script << "	legend.txt<-c(paste(\"GC\",sprintf (\"%.2f%%\",gcAvg),\"±\",sprintf (\"%.2f\",gcStd)))\n";
	script << "	legend('topright',legend.txt,bty='n')\n";
	script << "	if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)\n";
	script << "        }\n";
	script << "\n";
	script << "	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0.75,1), mar=c(3, 2, 2, 2),new=TRUE,cex.main=1)\n";
	script << "	legend.txt<-c(paste(\"A\",sprintf (\"%.2f%%\",aAvg),\"±\",sprintf (\"%.2f\",aStd)))\n";
	script << "	plot(Aaggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab=\"\",ylab=\"\",)\n";
	script << "\n";
	script << "	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0.5,0.75),mar=c(3, 2, 2, 2),new=TRUE)\n";
	script << "	legend.txt<-c(paste(\"T\",sprintf (\"%.2f%%\",tAvg),\"±\",sprintf (\"%.2f\",tStd)))\n";
	script << "	plot(Taggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab=\"\",ylab=\"\")\n";
	script << "\n";
	script << "	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0.25,0.5),mar=c(3, 2, 2, 2),new=TRUE)\n";
	script << "	legend.txt<-c(paste(\"C\",sprintf (\"%.2f%%\",cAvg),\"±\",sprintf (\"%.2f\",cStd)))\n";
	script << "	plot(Caggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab=\"\",ylab=\"\")\n";
	script << "\n";
	script << "	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0,0.25),mar=c(3, 2, 2, 2),new=TRUE)\n";
	script << "	legend.txt<-c(paste(\"G\",sprintf (\"%.2f%%\",gAvg),\"±\",sprintf (\"%.2f\",gStd)))\n";
	script << "	plot(Gaggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab=\"\",ylab=\"\")\n";
	script << "}\n";
	script << "if(file.exists(\"" << pre_nuc_composition_file << "\")){\n";
	script << "    readGC_plot(\"" << pre_nuc_composition_file << "\",qa.readsCount,0,0.5,\"Input Reads\",\"Number of reads (millions)\",FALSE)\n";
	script << "    readGC_plot(\"" << post_nuc_composition_file << "\",readsCount,0.5,1,\"Trimmed Reads\",\"\",TRUE)\n";
	script << "    #abline(v=0.5,lty=2,xpd=TRUE)\n";
	script << "} else {\n";
	script << "    readGC_plot(\"" << post_nuc_composition_file << "\",readsCount,0,1,\"\",\"Number of reads (millions)\",FALSE)\n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title(\"Reads GC content\",adj=0)\n";
	script << "\n";
	script << "\n";
	script << "#ATCG composition per base ATCG plot\n";
	script << "baseM<-read.table(file=\"" << post_base_matrix_file << "\")\n";
	script << "ATCG_composition_plot <- function(base_matrix_file,totalReads,xlab,ylab) {\n";
	script << "	baseM<-read.table(file=base_matrix_file)\n";
	script << "	aBase<-baseM$V1\n";
	script << "	tBase<-baseM$V2\n";
	script << "	cBase<-baseM$V3\n";
	script << "	gBase<-baseM$V4\n";
	script << "	nBase<-baseM$V5\n";
	script << "\n";
	script << "     rowSum<-rowSums(baseM)\n";
	script << "     rowSum[rowSum==0]<-1\n";
	script << "     aPer<-(aBase/rowSum)*100\n";
	script << "     tPer<-(tBase/rowSum)*100\n";
	script << "     cPer<-(cBase/rowSum)*100\n";
	script << "     gPer<-(gBase/rowSum)*100\n";
	script << "\n";
	script << "     ymax<-floor(max(aPer,tPer,cPer,gPer))\n";
	script << "     ymin<-floor(min(aPer,tPer,cPer,gPer))\n";
	script << "     if((ymin-5)>0){ymin <- ymin-5}else{ymin<-0}\n";
	script << "	xpos<-seq(1,length(aBase),1)\n";
	script << "	plot(xpos,aPer,col='green3',type='l',xaxt='n',xlab=xlab,ylab=ylab ,ylim=c(ymin,ymax+5))\n";
	script << "	lines(xpos,tPer,col='red')\n";
	script << "	lines(xpos,cPer,col='blue')\n";
	script << "	lines(xpos,gPer,col='black')\n";
	script << "	axis(1,at=xpos,labels=xpos)\n";
	script << "	legend('topright',c('A','T','C','G'),col=c('green3','red','blue','black'),box.col=0,lwd=1)\n";
	script << "	if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)\n";
	script << "        }\n";
	script << "	return(nBase)\n";
	script << "}\n";
	script << "if(file.exists(\"" << pre_base_matrix_file << "\")){\n";
	script << "    par(mfrow=c(1,2),mar=c(5,6,4,2))\n";
	script << "    qa.nBase<-ATCG_composition_plot(\"" << pre_base_matrix_file << "\",qa.readsCount,\"Input Reads Base\",'Base content (%)')\n";
	script << "    nBase<-ATCG_composition_plot(\"" << post_base_matrix_file << "\",readsCount,\"Trimmed Reads Base\",\"\")\n";
	script << "}else{\n";
	script << "    nBase<-ATCG_composition_plot(\"" << post_base_matrix_file << "\",readsCount,\"\",'Base content (%)')\n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title(\"Nucleotide Content Per Cycle\")\n";
	script << "\n";
	script << "\n";
	script << "#N composition per Base plot\n";
	script << "N_composition_plot<-function(BaseArray,totalReads,xlab,ylab){\n";
	script << "  xpos<-seq(1,length(BaseArray),1)\n";
	script << "  plot(xpos,BaseArray/totalReads*1000000,col='red',type='l',xaxt='n',xlab=xlab,ylab=ylab,ylim=c(0,max(BaseArray/totalReads*1000000)))\n";
	script << "  axis(1,at=xpos,labels=xpos)\n";
	script << "  legend('topright',paste(\"Total bases: \",sum(BaseArray)),bty='n')\n";
	script << "  if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "      mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)\n";
	script << "  }\n";
	script << "}\n";
	script << "if (sum(nBase) > 0){\n";
	script << "  if(file.exists(\"" << pre_base_matrix_file << "\")){\n";
	script << "    par(mfrow=c(1,2),mar=c(5,6,4,2))\n";
	script << "    N_composition_plot(qa.nBase,qa.readsCount,\"Input reads Position\",\"N Base count per million reads\")\n";
	script << "    N_composition_plot(nBase,readsCount,\"Trimmed reads Position\",\"\")\n";
	script << "  }else{\n";
	script << "    N_composition_plot(nBase,readsCount,\"Position\",\"N Base count per million reads\")\n";
	script << "  } \n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title(\"N Nucleotide Content Per Cycle\")\n";
	script << "\n";
	script << "if(file.exists(\"" << kmer_rarefaction_file << "\")){\n";
	script << "par(mar=c(5,6,4,2))\n";
	script << "kmerfile<-read.table(file=\"" << kmer_rarefaction_file << "\")\n";
	script << "sampling_size<-sum(kmerfile$V1)\n";
	script << "sampling<-\"\"\n";
	script << "if(sampling_size< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "    sampling<-paste(\"(Sampling\",format(sampling_size/1000000,digit=3),\"M Reads)\")\n";
	script << "}\n";
	script << "cumSeqNum<-cumsum(kmerfile$V1);\n";
	script << "plot(cumSeqNum/1000000,kmerfile$V3/1000000,xlab=paste(\"Number of Sequence (million)\",sampling), ylab=\"Number of Distinct K-mer (million,k=" << m_opt.kmer << ")\",type=\'l\',lty=2)\n";
	script << "lines(cumSeqNum/1000000,kmerfile$V2/1000000,col='blue',lwd=2)\n";
	script << "title(\"Kmer Rarefaction Curve\")\n";
	script << "y<-kmerfile$V2/1000000\n";
	script << "x<-cumSeqNum/1000000\n";
	script << "lres<-lm(y~x)\n";
	script << "# y=ax+b\n";
	script << "a<-format(coef(lres)[[2]], digits = 2)\n";
	script << "b<-format(coef(lres)[[1]], digits = 2)\n";
	script << "par(def.par)#- reset to default\n";
	script << "\n";
	script << "}\n";
	script << "\n";
	script << "if(file.exists(\"" << kmer_histogram_file << "\")){\n";
	script << "par(mar=c(5,6,4,2))\n";
	script << "sampling<-\"\"\n";
	script << "if(sampling_size<" << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "    sampling<-paste(\"(Sampling\",format(sampling_size/1000000,digit=3),\"M Reads)\")\n";
	script << "}\n";
	script << "kmerHfile<-read.table(file=\"" << kmer_histogram_file << "\")\n";
	script << "barplot(kmerHfile$V2[3:length(kmerHfile$V1)],names.arg=kmerHfile$V1[3:length(kmerHfile$V1)],xlab=\"K-mer Count (k=" << m_opt.kmer << ")\",log=\'y\',ylab=\"Number of K-mer\",main=paste(\"K-mer Frequency Histogram\",sampling),col=\'darkcyan\',border=\'darkcyan\')\n";
	script << "total_kmer<-sum(kmerHfile$V2)\n";
	script << "legend('topleft',paste(\"Total: \",total_kmer),bty='n')\n";
	script << "par(fig=c(0.5,0.9,0.5,1), new=TRUE)\n";
	script << "barplot(kmerHfile$V2[3:length(kmerHfile$V1)],xlim=c(0,100),names.arg=kmerHfile$V1[3:length(kmerHfile$V1)],log=\'y\',col=\'darkcyan\',border=\'darkcyan\')\n";
	script << "par(fig=c(0,1,0,1))\n";
	script << "par(def.par)#- reset to default\n";
	script << "\n";
	script << "}\n";
	script << "\n";
	script << "# read avg quality count barplot \n";
	script << "quality_histogram<-function(qual_histogram_file,totalReads,xlab,ylab){\n";
	script << "	Qhist_file<-read.table(file=qual_histogram_file,header=TRUE)\n";
	script << "	cumulate<-cumsum(Qhist_file$readsNum)\n";
	script << "	par(mar=c(5,6,5,4))\n";
	script << "	if (missing(ylab)){ylab2<-\"\"} else {ylab2<-ylab}\n";
	script << "	plot(Qhist_file$Score,Qhist_file$readsNum/1000000,type='h',xlim=c(max(Qhist_file$Score),min(Qhist_file$Score)),xlab=xlab, ylab=ylab2,lwd=12,lend=2)\n";
	script << "	par(new=TRUE)\n";
	script << "	plot(Qhist_file$Score,cumulate/sum(Qhist_file$readsNum)*100,type='l',xlim=c(max(Qhist_file$Score),min(Qhist_file$Score)),yaxt='n',xaxt='n',ylab=\"\",xlab=\"\",col='blue',lwd=3)\n";
	script << "	axis(4,col='blue',col.ticks='blue',col.axis='blue')\n";
	script << "	if (missing(ylab)){\n";
	script << "	  mtext(side=4,'Cumulative Percentage',line=2,col='blue')\n";
	script << "	}\n";
	script << "	Qover20Reads<-sum(as.numeric(Qhist_file$readsNum[Qhist_file$Score>=20]))\n";
	script << "	Qover20ReadsPer<-sprintf(\"%.2f%%\",Qover20Reads/sum(Qhist_file$readsNum)*100)\n";
	script << "	Qover20Bases<-sum(as.numeric(Qhist_file$readsBases[Qhist_file$Score>=20]))\n";
	script << "	Qover20AvgLen<-sprintf(\"%.2f\",Qover20Bases/Qover20Reads)\n";
	script << "	TotalBases<-sum(as.numeric(Qhist_file$readsBases))\n";
	script << "	mtext(side=3,paste(\"Number of Q>=20 reads:\",formatC(Qover20Reads,format=\"d\",big.mark=\",\"),\"(\",Qover20ReadsPer,\")\",\", mean Length:\",Qover20AvgLen),adj=0,cex=0.8,line=0.3)\n";
	script << "	if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "           mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),line=1,adj=0,cex=0.8)\n";
	script << "        }\n";
	script << "    return(TotalBases)\n";
	script << "}\n";
	script << "if(file.exists(\"" << pre_quality_histogram_file << "\"))\n";
	script << "{\n";
	script << "    par(mfrow=c(1,2))\n";
	script << "    qa.totalBases<-quality_histogram(\"" << pre_quality_histogram_file << "\",qa.readsCount,\"Input Reads Avg Score\",\"Reads Number (millions)\")\n";
	script << "    totalBases<-quality_histogram(\"" << post_quality_histogram_file << "\",readsCount,\"Trimmed Reads Avg Score\")\n";
	script << "}else{\n";
	script << "    totalBases<-quality_histogram(\"" << post_quality_histogram_file << "\",readsCount,\"Avg Score\",\"Reads Number (millions)\")\n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title('Reads Average Quality Histogram')\n";
	script << "\n";
	script << "# read in matrix file for the following three plots\n";
	script << "quality_boxplot<-function(quality_matrix_file,totalReads,totalBases,xlab,ylab){\n";
	script << "	z<-as.matrix(read.table(file=quality_matrix_file))\n";
	script << "	x<-1:nrow(z)\n";
	script << "	y<-1:ncol(z)\n";
	script << "	y<-y-1\n";
	script << "\n";
	script << "	#quality boxplot per base\n";
	script << "	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol\n";
	script << "	plot(1:length(x),x,type='n',xlab=xlab,ylab=ylab, ylim=c(0,max(y)+1),xaxt='n')\n";
	script << "	axis(1,at=x,labels=x)\n";
	script << "\n";
	script << "	for (i in 1:length(x)) {\n";
	script << "	  total<-sum(z[i,])\n";
	script << "	  qAvg<-sum(y*z[i,])/total\n";
	script << "	  if (is.wholenumber(total/2))\n";
	script << "	  {\n";
	script << "		 med<-( min(y[cumsum((z[i,]))>=total/2]) + min(y[cumsum((z[i,]))>=total/2+1]) )/2\n";
	script << "	  }\n";
	script << "	  else\n";
	script << "	  {\n";
	script << "		 med<-min(y[cumsum((z[i,]))>=ceiling(total/2)])\n";
	script << "	  }\n";
	script << "\n";
	script << "	  if (is.wholenumber(total/4))\n";
	script << "	  {\n";
	script << "		 Q1<-( min(y[cumsum((z[i,]))>=total/4]) + min(y[cumsum((z[i,]))>=total/4+1]) )/2\n";
	script << "	  }\n";
	script << "	  else\n";
	script << "	  {\n";
	script << "		 Q1<-min(y[cumsum((z[i,]))>=round(total/4)])\n";
	script << "	  }\n";
	script << "\n";
	script << "	  if (is.wholenumber(total/4*3))\n";
	script << "	  {\n";
	script << "		 Q3<-( min(y[cumsum((z[i,]))>=total/4*3]) + min(y[cumsum((z[i,]))>=total/4*3+1]) )/2\n";
	script << "	  }\n";
	script << "	  else\n";
	script << "	  {\n";
	script << "		 Q3<-min(y[cumsum((z[i,]))>=round(total/4*3)])\n";
	script << "	  }\n";
	script << "	  maxi<-max(y[z[i,]>0])\n";
	script << "	  mini<-min(y[z[i,]>0])\n";
	script << "	  #if (Q1 == 'Inf') {Q1 = maxi}\n";
	script << "	  if (Q3 == 'Inf') {Q3 = maxi}\n";
	script << "	  IntQ<-Q3-Q1\n";
	script << "	  mini<-max(mini,Q1-1.5*IntQ)\n";
	script << "	  maxi<-min(maxi,Q3+1.5*IntQ)\n";
	script << "	  rect(i-0.4,Q1,i+0.4,Q3,col='bisque')\n";
	script << "	  lines(c(i,i),c(Q3,maxi),lty=2)\n";
	script << "	  lines(c(i,i),c(mini,Q1),lty=2)\n";
	script << "	  lines(c(i-0.4,i+0.4),c(mini,mini))\n";
	script << "	  lines(c(i-0.4,i+0.4),c(maxi,maxi))\n";
	script << "	  lines(c(i-0.4,i+0.4),c(med,med))\n";
	script << "	  #points(i,qAvg,col='red')\n";
	script << "	  reads_num<-formatC(totalReads,format=\"d\",big.mark=\",\")\n";
	script << "	  reads_base<-formatC(totalBases,format=\"d\",big.mark=\",\")\n";
	script << "	  abline(h=20, col = \"gray60\")\n";
	script << "	  legend(\"bottomleft\",c(paste(\"# Reads: \",reads_num),paste(\"# Bases:\",reads_base)))\n";
	script << "	  if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "              mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)\n";
	script << "          }\n";
	script << "	## for outliers\n";
	script << "	#points()\n";
	script << "	}\n";
	script << "}\n";
	script << "if (file.exists(\"" << pre_quality_matrix_file << "\")){\n";
	script << "    par(mfrow=c(1,2),mar=c(5,6,4,2))\n";
	script << "    quality_boxplot(\"" << pre_quality_matrix_file << "\",qa.readsCount,qa.totalBases,\"Input Reads Position\",\"Quality score\")\n";
	script << "    quality_boxplot(\"" << post_quality_matrix_file << "\",readsCount,totalBases,\"Trimmed Reads Position\",\"\")\n";
	script << "}else{\n";
	script << "    quality_boxplot(\"" << post_quality_matrix_file << "\",readsCount,totalBases,\"Position\",\"Quality score\")\n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title(\"Quality Boxplot Per Cycle\")\n";
	script << "\n";
	script << "#quality 3D plot\n";
	script << "quality_3d_plot<-function(quality_matrix_file,totalReads,xlab,ylab){\n";
	script << "	z<-as.matrix(read.table(file=quality_matrix_file))\n";
	script << "	x<-1:nrow(z)\n";
	script << "	y<-1:ncol(z)\n";
	script << "	y<-y-1\n";
	script << "    persp(x,y,z/1000000,theta = 50, phi = 30, expand = 0.7, col = rev(terrain.colors(length(z),alpha=0.8)),border=NA,ntick=10,ticktype=\"detailed\",xlab=xlab,ylab=ylab,zlab=\"\",r=6,shade=0.75)\n";
	script << "    mtext(side=2, \"Frequency (millions)\",line=2)\n";
	script << "	if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)\n";
	script << "        }\n";
	script << "}\n";
	script << "if (file.exists(\"" << pre_quality_matrix_file << "\")){\n";
	script << "    par(mfrow=c(1,2),mar=c(5,6,4,2))\n";
	script << "    quality_3d_plot(\"" << pre_quality_matrix_file << "\",qa.readsCount,\"Input Reads Position\",\"Q Score\")\n";
	script << "    quality_3d_plot(\"" << post_quality_matrix_file << "\",readsCount,\"Trimmed Reads Position\",\"Q Score\")\n";
	script << "}else{\n";
	script << "    quality_3d_plot(\"" << post_quality_matrix_file << "\",readsCount,\"Position\",\"Q Score\")\n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title(\"Quality 3D plot. (Position vs. Score vs. Frequency)\")\n";
	script << "\n";
	script << "#Quality count bar plot\n";
	script << "upper_limit<-41\n";
	script << "quality_count_histogram<-function(quality_matrix_file,totalReads,highestScore,xlab,ylab){\n";
	script << "    z<-as.matrix(read.table(file=quality_matrix_file));\n";
	script << "    col<-colSums(z)\n";
	script << "    less30columnNum<-length(col)-highestScore+30-1\n";
	script << "    atleast30columnNum<-highestScore-30+1\n";
	script << "    color<-c(rep('blue',less30columnNum),rep('darkgreen',atleast30columnNum))\n";
	script << "    over30per<-sprintf(\"%.2f%%\",sum(col[(less30columnNum+1):length(col)])/sum(col)*100)\n";
	script << "    countInM<-col/1000000\n";
	script << "    avgQ<-sprintf(\"%.2f\",sum(seq(0,41,1)*col)/sum(col))\n";
	script << "    plot(seq(0,highestScore,1),countInM,col=color,type='h',ylab=ylab,xlab=xlab,lwd=12,lend=2,bty='n')\n";
	script << "    abline(v=29.5,col='darkgreen')\n";
	script << "    text(30,(max(countInM)-min(countInM))*0.9,labels=\">=Q30\",cex=0.8,adj=0,col='darkgreen')\n";
	script << "    text(30,(max(countInM)-min(countInM))*0.85,labels=over30per,cex=0.8,adj=0,col='darkgreen')\n";
	script << "    mtext(side=3,paste(\"Average: \",avgQ),adj=0)\n";
	script << "	if (totalReads< " << m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER] << "){\n";
	script << "            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),line=1 ,adj=0)\n";
	script << "        }\n";
	script << "}\n";
	script << "if (file.exists(\"" << pre_quality_matrix_file << "\")){\n";
	script << "    par(mfrow=c(1,2),mar=c(5,6,4,2))\n";
	script << "    quality_count_histogram(\"" << pre_quality_matrix_file << "\",qa.readsCount,upper_limit,\"Input Reads Q score\",\"Total (million)\")\n";
	script << "    quality_count_histogram(\"" << post_quality_matrix_file << "\",readsCount,upper_limit,\"Trimmed Reads Q score\",\"\")\n";
	script << "}else{\n";
	script << "    quality_count_histogram(\"" << post_quality_matrix_file << "\",readsCount,upper_limit,\"Q score\",\"Total (million)\")\n";
	script << "}\n";
	script << "par(def.par)#- reset to default\n";
	script << "title(\"Quality report\")\n";
	script << "\n";
	script << "tmp<-dev.off()\n";
	script << "\n";
	script << "quit()\n";
	
	// Instead of writing a temp file to save the R script, directly pipe the R script
	// to R. This avoids a potential security issue that might arise if the temporary file
	// was modified on disk after writing, but before executing.
	FILE *fpipe = popen("R --vanilla --silent --slave", "w");

	if(fpipe == NULL){
		cerr << "Warning: Unable to run R for plot generation" << endl;
	}
	else{
		fwrite(script.str().c_str(), sizeof(char), script.str().size(), fpipe);
		pclose(fpipe);
	}
	
	if(!m_opt.debug){
		// Clean up
		unlink( pre_quality_matrix_file.c_str() );
		unlink( post_quality_matrix_file.c_str() );

		unlink( pre_base_matrix_file.c_str() );
		unlink( post_base_matrix_file.c_str() );

		unlink( pre_quality_histogram_file.c_str() );
		unlink( post_quality_histogram_file.c_str() );

		unlink( pre_nuc_composition_file.c_str() );
		unlink( post_nuc_composition_file.c_str() );

		unlink( pre_length_histogram_file.c_str() );
		unlink( post_length_histogram_file.c_str() );

		unlink( kmer_rarefaction_file.c_str() );
		//unlink( kmer_files.c_str() );
		unlink( kmer_histogram_file.c_str() );
	}
}

void write_base_content(const string &m_filename, const vector<NucleotideCount> &m_data)
{
	ofstream fout( m_filename.c_str() );
	
	if(!fout){
	
		cerr << "Warning: Unable to open " << m_filename << " for writing base content data!" << endl;
		return;
	}
	
	fout << setprecision(2) << fixed;
	
	// Composition for nucleotide A:
	for(unsigned int i = 0;i < NUM_COMPOSITION_BIN;++i){
		
		const size_t &num = m_data[i].num_A;
		
		if(num != 0){
			fout << "A\t" << i*0.01 << '\t' << num << endl;
		}
	}
	
	// Composition for nucleotide T:
	for(unsigned int i = 0;i < NUM_COMPOSITION_BIN;++i){
		
		const size_t &num = m_data[i].num_T;
		
		if(num != 0){
			fout << "T\t" << i*0.01 << '\t' << num << endl;
		}
	}
	
	// Composition for nucleotide C:
	for(unsigned int i = 0;i < NUM_COMPOSITION_BIN;++i){
		
		const size_t &num = m_data[i].num_C;
		
		if(num != 0){
			fout << "C\t" << i*0.01 << '\t' << num << endl;
		}
	}
	
	// Composition for nucleotide G:
	for(unsigned int i = 0;i < NUM_COMPOSITION_BIN;++i){
		
		const size_t &num = m_data[i].num_G;
		
		if(num != 0){
			fout << "G\t" << i*0.01 << '\t' << num << endl;
		}
	}
	
	// Composition for nucleotide N:
	for(unsigned int i = 0;i < NUM_COMPOSITION_BIN;++i){
		
		const size_t &num = m_data[i].num_N;
		
		if(num != 0){
			fout << "N\t" << i*0.01 << '\t' << num << endl;
		}
	}
	
	// Composition for nucleotide GC:
	for(unsigned int i = 0;i < NUM_COMPOSITION_BIN;++i){
		
		const size_t &num = m_data[i].num_GC;
		
		if(num != 0){
			fout << "GC\t" << i*0.01 << '\t' << num << endl;
		}
	}
}

void write_matrix(const string &m_filename, const matrix<size_t> &m_data)
{
	const size_t num_row = m_data.get_num_row();
	const size_t num_col = m_data.get_num_col();
	
	if( (num_row == 0) || (num_col == 0) ){
		return;
	}
	
	ofstream fout( m_filename.c_str() );
	
	if(!fout){
	
		cerr << "Warning: Unable to open " << m_filename << " for writing matrix data!" << endl;
		return;
	}
	
	for(size_t i = 0;i < num_row;++i){
		
		fout << m_data(i, 0);
		
		for(size_t j = 1;j < num_col;++j){
			fout << '\t' << m_data(i, j);
		}
		
		fout << endl;
	}
}

void write_quality_histogram(const string &m_filename, const vector<size_t> &m_read_hist,
	const vector<size_t> &m_base_hist)
{
	ofstream fout( m_filename.c_str() );
	
	if(!fout){
	
		cerr << "Warning: Unable to write quality histogram file: " << m_filename << endl;
		return;
	}
	
	if( (m_read_hist.size() <= MAX_QUALITY_SCORE) || 
	    (m_base_hist.size() <= MAX_QUALITY_SCORE) ){
		throw __FILE__ ":write_quality_histogram: Histogram is too small!";
	}
	
	fout << "Score\treadsNum\treadsBases" << endl;
	
	for(int i = MAX_QUALITY_SCORE;i >= 0;--i){
		fout << i << '\t' << m_read_hist[i] << '\t' << m_base_hist[i] << endl;
	}
}

void write_length_histogram(const string &m_filename, const vector<size_t> &m_hist)
{
	ofstream fout( m_filename.c_str() );
	
	if(!fout){
	
		cerr << "Warning: Unable to write length histogram file: " << m_filename << endl;
		return;
	}
	
	const size_t len = m_hist.size();
	
	// Start at a length of 1 (not zero)
	for(size_t i = 1;i < len;++i){
		fout << i << '\t' << m_hist[i] << endl;
	}
}

void write_kmer_histogram(const string &m_filename, const MAP<size_t, size_t> &m_hist)
{
	ofstream fout( m_filename.c_str() );

        if(!fout){

                cerr << "Warning: Unable to write kmer histogram file: " << m_filename << endl;
                return;
        }

	// Output the kmer frequency histogram in order of ascending frequency value
	vector<size_t> f;

	f.reserve( m_hist.size() );

	for(MAP<size_t, size_t>::const_iterator i = m_hist.begin();i != m_hist.end();++i){
		f.push_back(i->first);
	}

	sort( f.begin(), f.end() );

	for(vector<size_t>::const_iterator i = f.begin();i != f.end();++i){
		
		MAP<size_t, size_t>::const_iterator iter = m_hist.find(*i);

		if( iter == m_hist.end() ){
			throw __FILE__ ":write_kmer_histogram: Unable to look up kmer frequency";
		}

		fout << *i << ' ' << iter->second << endl;
	}
}

void write_rarefaction(const string &m_filename, const vector<Rarefaction> &m_rarefaction)
{
	ofstream fout( m_filename.c_str() );

        if(!fout){

                cerr << "Warning: Unable to write kmer rarefaction file: " << m_filename << endl;
                return;
        }

	size_t last_num_seq = 0;

	for(vector<Rarefaction>::const_iterator i = m_rarefaction.begin();i != m_rarefaction.end();++i){

		fout << i->num_seq - last_num_seq << '\t' << i->distinct_kmer << '\t' << i->total_kmer << endl;	
		last_num_seq = i->num_seq;
	}
}
