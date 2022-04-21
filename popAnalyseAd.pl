#! /usr/bin/perl -w
#
use strict;
use Getopt::Long;
use Term::ANSIColor;
use File::Basename;


my ($pheno,$pop,$vcf,$pname,$odir,$setNum,$pcanum,$maf,$geno,$mind,$hwe,$mac,$minQ,$minDP);
my $help;

GetOptions(
        "pt:s"            =>\$pheno,
        "pop:s"           =>\$pop,
        "vcf:s"           =>\$vcf,
        "pn:s"            =>\$pname,

        "n:i"             =>\$setNum,
        "pca:i"           =>\$pcanum,
        "maf:s"           =>\$maf,
        "geno:s"          =>\$geno,
        "mind:s"          =>\$mind,
        "hwe:s"           =>\$hwe,

        "mac:i"           =>\$mac,
        "minQ:i"          =>\$minQ,
        "minDP:i"          =>\$minDP,

        "o:s"            =>\$odir,

        'h'               =>\$help,
);

sub usage{
        my $perl_name=basename $0;
        print color 'bold blue';
        print "
--------------------------------------------------------------------------------------------------------
        This process is used for population analysis, including Kship,PCA,Tree,Gwas and structure

Usage:
                -pt             pheno:sampleName\\tsampleName\\tphenotype value
                -pop            population information:sampleName\\tsub population\\tsub population
                -vcf            vcf format file
                -pn             phenotype name:test
                -n              chr num:20
                -pca            PCA:3
                -maf            maf:0.05
                -geno           geno:0.2
                -mind           mind:0.2
                -hwe            hwe:NA
                -mac            mac:3
                -minQ           minQ:20
                -minDP          minDP:3
                -o              out dir

                -h              help

                 eg:
                        $perl_name

Note:

        If you have any questions ,please contact to super_qd_shao\@163.com
--------------------------------------------------------------------------------------------------------
                                Good Good Study,Day Day Up
........................................................................................................
\n";
        exit(0);
}
if ($help){&usage();}
unless(defined $odir && defined $vcf ){&usage();}
mkdir $odir unless -d $odir;
$pname  ||="test";
$setNum ||=20;
$pcanum ||=3;
$maf    ||=0.05;
$geno   ||=0.2;
$mind   ||=0.2;
$hwe    ||='';
$mac    ||=3;
$minQ   ||=20;
$minDP  ||=3;

my $part=basename $pheno;
print "\n\n#----------------------------------PART:$part-------------------------------------------#\n";
print "#----------------------------------------------------------------------------------#\n";
print "pheno : $pheno\n";
print "pop   : $pop\n";
print "vcf   ：$vcf\n";
print "pname : $pname\n";
print "odir  ：$odir\n";
print "setNum：$setNum\n";
print "pcanum：$pcanum\n";
print "#----------------------------------------------------------------------------------#\n";
print "#-------------------Everyone is always gonna go for the better everything.---------#\n";

my $Bin="./BGISEQ_WGRS/";

my $gwasdir="$odir/05.Gwas";mkdir $gwasdir unless -d $gwasdir;
my $workdir="$gwasdir/$pname";mkdir $workdir unless -d $workdir;
my $datadir="$odir/01.data";mkdir $datadir unless -d $datadir;
my $beddir="$workdir/r0.bed";mkdir $beddir unless -d $beddir;
my $qcdir="$odir/00.QC";mkdir $qcdir unless -d $qcdir;

open  QC,">$qcdir/QC.sh" or die $!;
print QC "sh $Bin/bin/QC/QC.sh $vcf $pname $qcdir\n";
close QC;

open   STEP0,">$datadir/work.filter.sh" or die $!;
print  STEP0 "$Bin/bin/vcftools --gzvcf  $vcf --max-missing 0.9  --maf 0.1 --max-maf 0.5 --mac $mac --minQ $minQ --minDP $minDP  --recode --out $datadir/result\n";
if($hwe){
        print  STEP0 "$Bin/bin/plink --allow-extra-chr --chr-set $setNum --hwe $hwe --memory 5000  --recode  --make-bed  --out $datadir/result --vcf $datadir/result.recode.vcf\n";
}else{
        print  STEP0 "$Bin/bin/plink --allow-extra-chr --chr-set $setNum --memory 5000  --recode  --make-bed  --out $datadir/result --vcf $datadir/result.recode.vcf\n";
}

print  STEP0 "mv $datadir/result.fam $datadir/result.fam.tmp\nperl $Bin/bin/checkSampleName.pl $datadir/result.fam.tmp >$datadir/result.fam\n";
print  STEP0  "mv $datadir/result.map $datadir/result.1.map \nmv $datadir/result.bim $datadir/result.1.bim\n";

print  STEP0 "perl -ne 'chomp;my \@F=split;\$F[1]=\"\$F[0]:\$F[3]\";my \$out=join \"\\t\",\@F;print \"\$out\\n\";' $datadir/result.1.map >$datadir/result.map\n";
print  STEP0 "perl -ne 'chomp;my \@F=split;\$F[1]=\"\$F[0]:\$F[3]\";my \$out=join \"\\t\",\@F;print \"\$out\\n\";' $datadir/result.1.bim >$datadir/result.bim\n";

print  STEP0 "$Bin/bin/plink --allow-extra-chr --chr-set $setNum --memory 5000  --make-bed  --bfile $datadir/result  --recode 12  --transpose --output-missing-genotype 0 --out $datadir/result\n";
close  STEP0;

open   STEP1,">$beddir/work.tped.sh" or die $!;
print  STEP1 "$Bin/bin/plink --bfile $datadir/result --chr-set $setNum  --memory 5000  --geno $geno  --maf $maf   --mind $mind --biallelic-only  --make-bed --out $beddir/$pname --recode --allow-extra-chr\n";
print  STEP1  "mv $beddir/$pname.map $beddir/$pname.1.map \nmv $beddir/$pname.bim $beddir/$pname.1.bim\n";
print  STEP1 "perl -ne 'chomp;my \@F=split;\$F[1]=\"\$F[0]:\$F[3]\";my \$out=join \"\\t\",\@F;print \"\$out\\n\";' $beddir/$pname.1.map >$beddir/$pname.map\n";
print  STEP1 "perl -ne 'chomp;my \@F=split;\$F[1]=\"\$F[0]:\$F[3]\";my \$out=join \"\\t\",\@F;print \"\$out\\n\";' $beddir/$pname.1.bim >$beddir/$pname.bim\n";
print  STEP1 "$Bin/bin/plink --bfile $beddir/$pname   --chr-set $setNum --memory 5000  --pheno $pheno --allow-extra-chr --recode 12  --transpose --output-missing-genotype 0 --out $beddir/$pname\n";
print  STEP1 "$Bin/bin/plink --tfile $beddir/$pname    --chr-set $setNum --memory 5000  --pheno $pheno --recode  --allow-extra-chr --make-bed  --out $beddir/$pname\n";
print  STEP1 "mv $beddir/$pname.tfam $beddir/$pname.tfam.tmp\nperl $Bin/bin/checkSampleName.pl $beddir/$pname.tfam.tmp >$beddir/$pname.tfam\n";
close  STEP1;



kinship();
pca    ();
tree   ();
structure ();
gwas   ();
print "#----------------------------------------------------------------------------------#\n";

sub kinship{
        my $kinshipdir="$odir/02.Kship";
        mkdir $kinshipdir unless -d $kinshipdir;

        open  KINSHIP,">$kinshipdir/work.kinship.sh" or die $!;
        print KINSHIP "$Bin/bin/emmax-kin-intel64 -v  -s -d 10 $datadir/result -o $kinshipdir/kinship.hIBS.kinf\n";
        print KINSHIP "$Bin/bin/emmax-kin-intel64 -v  -d 10 $datadir/result  -o $kinshipdir/kinship.hBN.kinf\n";
        print KINSHIP "Rscript $Bin/bin/kinship.hIBS.R $kinshipdir\n";
        print KINSHIP "Rscript $Bin/bin/kinship.hBN.R   $kinshipdir\n";
        close KINSHIP;
}



sub pca{
        my $pcadir="$odir/03.PCA";
        my $gcta="$pcadir/r1.gcta";
        my $smartpca="$pcadir/r2.EIG";
        mkdir $pcadir unless -d $pcadir;
        mkdir $gcta   unless -d  $gcta;
        mkdir $smartpca unless -d $smartpca;
        my $pca_pop=$pop;

        open  PCA1,">$gcta/work.gcta.sh" or die $!;
        print PCA1 "$Bin/bin/plink  --tfile $datadir/result --chr-set $setNum --pheno $pheno  --pca $pcanum --out $gcta/pca.plink\n";
        print PCA1 "$Bin/bin/gcta64 --bfile $datadir/result  --make-grm --pca $pcanum --autosome-num  $setNum --autosome --out $gcta/gcta.plink.result_grm\n";
        print PCA1 "$Bin/bin/gcta64 --grm $gcta/gcta.plink.result_grm   --pca $pcanum --autosome-num  $setNum  --out $gcta/gcta.plink.result_pca\n";
        print PCA1 "$Bin/bin/Rscript $Bin/bin/PCA.R  $gcta/gcta.plink.result_pca   $pca_pop  $gcta/gcta_pca.result\n";
        close PCA1;


        open  PCA2,">$smartpca/work.smartpca.sh" or die $!;
        print PCA2 "ln -s $datadir/result* $smartpca\n";
        print PCA2 "##convert plink to EIG\n#create .conf file first\n";
        print PCA2 "echo \"typename:         result.ped\" > $smartpca/transfer.conf\n";
        print PCA2 "echo \"snpname:          result.map\" >> $smartpca/transfer.conf\n";
        print PCA2 "echo \"indivname:        result.ped\" >> $smartpca/transfer.conf\n";
        print PCA2 "echo \"outputformat:     EIGENSTRAT\" >> $smartpca/transfer.conf\n";
        print PCA2 "echo \"genotypeoutname:  smartpca.geno\" >> $smartpca/transfer.conf\n";
        print PCA2 "echo \"snpoutname:       smartpca.snp\" >> $smartpca/transfer.conf\n";
        print PCA2 "echo \"indivoutname:     smartpca.ind\" >> $smartpca/transfer.conf\n";
        print PCA2 "$Bin/bin/convertf -p $smartpca/transfer.conf >$smartpca/convert.log &&\n";
        print PCA2 "awk \'\{print \$1,\$2,\"Control\"\}\' $smartpca/smartpca.ind >$smartpca/smartpca_update.ind &&\n";
        print PCA2 "mv $smartpca/smartpca_update.ind $smartpca/smartpca.ind&&\n";
        print PCA2 "##convert EIG to smartpca\n#create .conf file first\n";
        print PCA2 "echo \"genotypename:     smartpca.geno\" > $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"snpname:          smartpca.snp\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"indivname:        smartpca.ind\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"evecoutname:      smartpca.eigenvec\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"evaloutname:      smartpca.eigenval\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"altnormstyle:     NO\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"numoutevec:       10\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"numoutlieriter:   5\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "echo \"outliersigmathresh: 6.0\" >> $smartpca/runningPCA.conf\n";
        print PCA2 "$Bin/bin/smartpca -p $smartpca/runningPCA.conf >$smartpca/runningPCA.log &&\n";
        print PCA2 "sed -i \'2,\$s/:/\t/\' $smartpca/smartpca.eigenvec &&\n";
        print PCA2 "$Bin/bin/Rscript $Bin/bin/PCA.R  $smartpca/smartpca $pca_pop  $smartpca/smartpca.result\n";
        close PCA2;
}


sub tree{
        my $treedir="$odir/04.Tree";
        my $SNPhylo="$treedir/r1.SNPhylo";
        my $treebest="$treedir/r2.treebest";
        mkdir $treedir unless -d $treedir;
        mkdir $SNPhylo unless -d $SNPhylo;
        mkdir $treebest unless -d $treebest;
        my $tree_pop=$pop;

        my $fasta="$datadir/result.recode.vcf";$fasta=~s/\.vcf*/\.fasta/;
        open  TREE1,">$SNPhylo/work.ml.sh" or die $!;
        print TREE1 "$Bin/bin/jre1.8.0_101/bin/java -Xmx1G  -XX:ParallelGCThreads=1 -jar $Bin/bin/NGSEPcore_3.3.0.jar ConvertVCF -printHapmap $datadir/result.recode.vcf  $SNPhylo/genotype\n";
        print TREE1 "sh $Bin/bin/SNPhylo/snphylo.sh -a $setNum  -H $SNPhylo/genotype_hmp.txt  -A -P  $SNPhylo/snphylo\n";
        print TREE1 "$Bin/bin/Rscript $Bin/bin/tree.R $SNPhylo/snphylo.ml.tree $tree_pop\n";
        close TREE1;

        open  TREE2,">$treebest/work.nj.sh" or die $!;
        print TREE2 "/dellfsqd2/ST_OCEAN/USER/hujie/soft/miniconda3/bin/python2 $Bin/bin/vcf2phylip.py -i $datadir/result.recode.vcf  -f -n -b\n";  ##combine.snp.min4.phy
        print TREE2 "$Bin/bin/treebest nj  $datadir/result.recode.fasta >$treebest/nj_tree.out\n";
        print TREE2 "$Bin/bin/Rscript $Bin/bin//tree.R $treebest/nj_tree.out $tree_pop\n";
        close TREE2;
}

sub gwas{
        my $plink  ="$workdir/r1.plink";
        mkdir   $plink   unless -d $plink;

        open  RGWAS,">$plink/gwas.R" or die $!;
        print RGWAS "library(CMplot)\n";
        print RGWAS "data=read.table(\"$plink/gwas.plink.qc.assco.result\")\n";
        print RGWAS "CMplot(data[,1:4],plot.type=\"b\",band=1,file=\"jpg\")\n";
        close RGWAS;
##PLINK
#The basic association commands (--assoc, --model, --fisher, --linear and --logistic) will test only a single phenotype. If your alternate phenotype file contains more than one phenotype, then adding the --all-pheno flag will make PLINK cycle over each phenotype, e.g. instead of a single plink.assoc output file, if there are 100 phenotypes, PLINK will now show
        open  GWAS,">$plink/work.plink.sh" or die $!;
        print GWAS "$Bin/bin/plink --tfile  $beddir/$pname --make-bed  --pheno $pheno   --chr-set $setNum  --out $plink/gwas.plink.qc\n";
        if(defined $pheno && -e $pheno ){
                #print GWAS "$Bin/bin/plink --bfile $plink/gwas.plink.qc --chr-set $setNum --pheno $pheno  --linear --allow-no-sex --out $plink/gwas.plink.qc.linear\n";
                print GWAS "$Bin/bin/plink --bfile $plink/gwas.plink.qc --chr-set $setNum --pheno $pheno   --allow-no-sex  --assoc --out $plink/gwas.plink.qc.assoc\n";
                print GWAS "$Bin/bin/plink --bfile $plink/gwas.plink.qc --chr-set $setNum --pheno $pheno  --covar $gwasdir/../p1.Kship/kinship.hIBS.kinf  --allow-no-sex  --assoc --adjust --gc --out  $plink/gwas.plink.qc.assoc.qc\n";
        #       print GWAS "$Bin/bin/plink --bfile $plink/gwas.plink.qc --chr-set $setNum --pheno $pheno  --covar $gwasdir/../p0.Kship/kinship.hIBS.kinf  --allow-no-sex  --logistic  --out $plink/gwplink.qc.logistic\n";


        }else{
                print GWAS "$Bin/bin/plink --bfile $plink/gwas.plink.qc --chr-set $setNum   -assoc --allow-no-sex --adjust --out $plink/gwas.plink.qc\n";

        }
        print GWAS "less $plink/gwas.plink.qc.assoc.qc.qassoc|awk '{print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$9}'|grep -v NA|sed -n 's/://p'  >$plink/gwas.plink.qc.assco.result\n";
        print GWAS "$Bin/bin/Rscript  $plink/gwas.R\n";
#       print GWAS "for file in *.pdf;do /usr/bin/convert \$file  \$file.png;done\n";
        close GWAS;

#EMMAX
#EMMAX is a statistical test for large scale human or model organism association mapping accounting for the sample structure. In addition to the computational efficiency obtained by EMMA algorithm, EMMAX takes advantage of the fact that each loci explains only a small fraction of complex traits, which allows us to avoid repetitive variance component estimation procedure, resulting in a significant amount of increase in computational time of association mapping using mixed model.
        my $emmax  ="$workdir/r2.emmax";
        mkdir   $emmax   unless -d $emmax;
        open  ER,">$emmax/gwas.R" or die $!;
        print ER "library(CMplot)\n";
        print ER "data=read.table(\"$emmax/gwas.result\")\n";
        print ER "CMplot(data)\n";
        print ER "CMplot(data,plot.type = \"m\",threshold = c(1e-6),threshold.col=c('grey'), threshold.lty = c(1),threshold.lwd = c(1), amplify = T,signal.cex = c(1), signal.pch = c(20),signal.col = c(\"red\"),file=\"jpg\",memo=\"1\",dpi=300)\n";##c(0.01,0.05)/nrow(data)
        print ER "CMplot(data,plot.type = \"m\",threshold = c(1e-6),threshold.col=c('grey'), threshold.lty = c(1),threshold.lwd = c(1), amplify = T,signal.cex = c(1), signal.pch = c(20),signal.col = c(\"red\"),file=\"pdf\",memo=\"1\",dpi=300)\n";
        ##print ER "CMplot(data,plot.type = \"m\",threshold = c(1e-6),threshold.col=c('grey','black'), threshold.lty = c(1,2),threshold.lwd = c(1,1), amplify = T,signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c(\"red\",\"orange\"),file=\"jpg\",memo=\"1\",dpi=300)\n";##c(0.01,0.05)/nrow(data)
        print ER "CMplot(data,plot.type = \"m\",amplify = T,signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c(\"red\",\"orange\"),file=\"jpg\",memo=\"2\",dpi=300)\n";

        print ER "CMplot(data,plot.type = \"d\",bin.size = 1e6, col = c(\"darkgreen\",\"yellow\",\"red\"))\n";
        print ER "CMplot(data,plot.type = \"q\",threshold = 0.05)\n";

        print ER "library(rMVP)\n";
        print ER "MVP.Report(data, plot.type=\"m\", LOG10=TRUE, ylim=NULL, threshold=c(1e-8,1e-6),threshold.lty=c(1,2), col=c(\"grey60\",\"grey30\"), threshold.lwd=c(1,1), threshold.col=c(\"black\",\"grey\"), amplify=TRUE, chr.den.col=c(\"darkgreen\", \"yellow\", \"red\"),bin.size=1e6,signal.col=c(\"red\",\"green\"),signal.cex=c(1,1),signal.pch=c(19,19),file.type=\"jpg\",memo=\"3\",dpi=300)\n";
##      print ER "MVP.Report(data, plot.type=\"m\", LOG10=TRUE, ylim=NULL,mplify=TRUE, chr.den.col=c(\"darkgreen\", \"yellow\", \"red\"),bin.size=1e6,signal.col=c(\"red\",\"green\"),signal.cex=c(1,1),signal.pch=c(19,19),file.type=\"jpg\",memo=\"4\",dpi=300)\n";

        #print ER "MVP.Report(data,plot.type=\"q\",conf.int.col=NULL,box=TRUE,file.type=\"jpg\",memo=\"\",dpi=300)\n";
        print ER "MVP.Data(fileBed=\"$beddir/$pname\",filePhe=\"$emmax/MVP.txt\",fileKin=TRUE,filePC=TRUE,out=\"$emmax/MVP\")\n";
        print ER "pheno <- read.table(\"$emmax/MVP.phe\", header = TRUE)\n";
        print ER "geno <- attach.big.matrix(\"$emmax/MVP.geno.desc\")\n";
        print ER "map <- read.table(\"$emmax/MVP.geno.map\", header = TRUE)\n";
        print ER "Kinship <- attach.big.matrix(\"$emmax/MVP.kin.desc\")\n";
        print ER "Covariates <- attach.big.matrix(\"$emmax/MVP.pc.desc\")\n";
        print ER "MVP.Hist(phe=pheno, file.type=\"jpg\", breakNum=18, dpi=300)\n";
        print ER "MVP.Report(map[, c(1:3)], plot.type=\"d\", col=c(\"darkgreen\", \"yellow\", \"red\"), file.type=\"jpg\", dpi=300)\n";
        print ER "pca <- attach.big.matrix(\"$emmax/MVP.pc.desc\")[, 1:3]\n";
        print ER "MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, col=c(\"red\", \"green\", \"yellow\"), file.type=\"jpg\")\n";
        close ER;
##
        open PCAR,">$emmax/pca.R" or die $!;
        print PCAR "library(\"ggplot2\")\n";
        print PCAR "a=read.table(\"$emmax/$pname.evec\",header=T)\n";
        print PCAR "ggplot(a,aes(PC1,PC2,color=species,pch = species))+geom_point(alpha=0.8,size=4)\n";
        print PCAR "ggplot(a,aes(PC1,PC3,color=species,pch = species))+geom_point(alpha=0.8,size=4)\n";
        print PCAR "ggplot(a,aes(PC2,PC3,color=species,pch = species))+geom_point(alpha=0.8,size=4)\n";
        close  PCAR;

        open  EMMAX,">$emmax/work.emmax.sh" or die $!;
        print EMMAX "cd $emmax\n";
        print EMMAX "less $pheno|awk '{print \$1\"\\t\"\$3}' >$emmax/MVP.txt\n";
        #print EMMAX "$Bin/bin/plink --bfile $beddir/$pname --pheno $pheno  --noweb  --make-bed --chr-set $setNum   -transpose --output-missing-genotype 0 --out $emmax/$pname\n";
        print EMMAX "$Bin/bin/emmax-kin-intel64 -v -d 10 $beddir/$pname -o $emmax/$pname.hBN.kinf -m 0.05 -c 0.05\n";
        #print EMMAX "perl  $Bin/bin/PCA.pl -i  $beddir/$pname.ped -a $beddir/$pname.map -b $pop   -o $emmax/$pname  -p $emmax/$pname.plot -e $emmax/$pname.evel -l $emmax/$pname.log\n";##PCA.evec
        print EMMAX "$Bin/bin/plink --bfile $beddir/$pname --maf 0.05 --geno 0.05  --chr-set $setNum  --allow-extra-chr --pca 10 --out $emmax/$pname\n";
        print EMMAX "awk 'BEGIN{a=1}{printf(\"\%s \%s \",\$1,\$1);printf(\"\%s \",a);for(i=3;i<12;i++){printf(\"\%s \",\$i)}printf(\"\%s\",\"\\n\")}' $emmax/$pname.eigenvec >$emmax/$pname.evec.PCA\n";
        ##eigenvec
        print EMMAX "perl $Bin/bin/SortSample.pl $beddir/$pname.nosex  $emmax/$pname.evec.PCA   >$emmax/$pname.evec.PCA.10\n";
        print EMMAX "perl $Bin/bin/SortSample.pl $beddir/$pname.nosex $pheno >$emmax/phenotype.txt\n";
###     print EMMAX "$Bin/bin/emmax -d 10 -t $beddir/$pname -p  $emmax/phenotype.txt  -c $emmax/$pname.evec.PCA.10  -k $emmax/$pname.hBN.kinf -o $emmax/$pname\n";
        print EMMAX "$Bin/bin/emmax -d 10 -t $beddir/$pname -p  $emmax/phenotype.txt   -k $emmax/$pname.hBN.kinf -o $emmax/$pname\n";
        print EMMAX "less $emmax/$pname.ps|sed -n 's/:/\\t/'p|awk '{print \$2\"\\t\"\$1\"\\t\"\$2\"\\t\"\$NF}' >$emmax/gwas.result\n";
        print EMMAX "Rscript $emmax/gwas.R\n";
#       print EMMAX "for file in *.pdf;do /usr/bin/convert \$file  \$file.png;done\n";
        close EMMAX;

        my $lmm    ="$workdir/r3.LMM";
        mkdir   $lmm     unless -d $lmm;
        open  LMM,">$lmm/work.lmm.sh" or die $!;
        print LMM "/zfsqd1/ST_OCEAN/USRS/shaolibin/tools/GWAS/FaSTLMM.207c.Linux/Linux_MKL/fastlmmc -tfile $beddir/$pname -tfileSim $beddir/$pname -pheno $pheno  -out $lmm/$pname.out.txt\n";
        print LMM "less  $lmm/$pname.out.txt|awk '{print \$1\"\\t\"\$2\"\\t\"\$4\"\\t\"\$5}'|grep -v NA|sed -n 's/://p'  >$lmm/gwas.plink.qc.assco.result\n";
        print LMM "Rscript  ./gwas.R\n";
        #print LMM "for file in *.pdf;do /usr/bin/convert \$file  \$file.png;done\n";
        close LMM;


        my $gemma    ="$workdir/r4.GEMMA";
        mkdir   $gemma unless -d $gemma;
        open GEMMA,">$gemma/work.gemma..sh" or die $!;
        open GEMMAR,">$gemma/gwas.R" or die $!;
        print GEMMAR "library(CMplot)\n";
        print GEMMAR "data=read.table(\"$gemma/gwas.GEMMA.result\")\n";
        print GEMMAR "CMplot(data[,1:4],plot.type=\"b\",band=1,file=\"jpg\",memo=\"GEMMA\",dpi=300)\n";
        close GEMMAR;

        print GEMMA "cd $gemma\n";
        print GEMMA "/zfsqd1/ST_OCEAN/USRS/shaolibin/tools/GWAS/Gemma/GEMMA/bin/gemma    -bfile $beddir/$pname -gk\n";
        print GEMMA "$Bin/bin/plink --tfile $beddir/$pname --chr-set $setNum  --allow-extra-chr --pca 10 --out $gemma/$pname\n";

        print GEMMA "awk 'BEGIN{a=1}{printf(\"\%s \",a);for(i=3;i<12;i++){printf(\"\%s \",\$i)}printf(\"\%s\",\"\\n\")}' $gemma/$pname.eigenvec >$gemma/$pname.evec.PCA\n";
        print GEMMA "perl $Bin/bin/SortSample.pl $pop  $gemma/$pname.evec.PCA   >$gemma/$pname.evec.PCA.10\n";

        print GEMMA "/zfsqd1/ST_OCEAN/USRS/shaolibin/tools/GWAS/Gemma/GEMMA/bin/gemma   -bfile $beddir/$pname  -lmm -k $gemma/output/result.cXX.txt -c $gemma/$pname.evec.PCA\n";
#       print GEMMA "/zfsqd1/ST_OCEAN/USRS/shaolibin/tools/GWAS/Gemma/GEMMA/bin/gemma   -bfile $beddir/$pname  -lmm -k $gemma/output/result.cXX.txt\n";
        print GEMMA "less  $gemma/output/result.assoc.txt|awk '{print \$2\"\\t\"\$1\"\\t\"\$3\"\\t\"\$NF}'|grep -v NA|sed -n 's/://p'  >$gemma/gwas.GEMMA.result\n";
        print GEMMA "Rscript $gemma/gwas.R\n";
        #print GEMMA "for file in *.pdf;do /usr/bin/convert \$file  \$file.png;done\n";
        close GEMMA;


        my $mvp    ="$workdir/r5.MVP";
        mkdir   $mvp unless -d $mvp;
#       `less $pheno|awk '{print \$1\"\\t\"\$3}' >$mvp/MVP.txt`;
        open GH,">$mvp/work.mvp.sh" or die $!;
        open GR,">$mvp/gwas.R" or die $!;
        print GR "library(rMVP)\n";
        print GR "MVP.Data(fileBed=\"$beddir/$pname\",filePhe=\"$mvp/MVP.txt\",fileKin=TRUE,filePC=TRUE,out=\"$mvp/mvp\")\n";
        print GR "pheno <- read.table(\"$mvp/mvp.phe\", header = TRUE)\n";
        print GR "geno <- attach.big.matrix(\"$mvp/mvp.geno.desc\")\n";
        print GR "map <- read.table(\"$mvp/mvp.geno.map\", header = TRUE)\n";
        print GR "Kinship <- attach.big.matrix(\"$mvp/mvp.kin.desc\")\n";
        print GR "Covariates <- attach.big.matrix(\"$mvp/mvp.pc.desc\")\n";

        print GR "MVP.Hist(phe=pheno, file.type=\"jpg\", breakNum=18, dpi=300)\n";
        print GR "MVP.Report(map[, c(1:3)], plot.type=\"d\", col=c(\"darkgreen\", \"yellow\", \"red\"), file.type=\"jpg\", dpi=300)\n";
        print GR "pca <- attach.big.matrix(\"$mvp/mvp.pc.desc\")[, 1:3]\n";
        print GR "MVP.PCAplot(PCA=pca, Ncluster=3, class=NULL, col=c(\"red\", \"green\", \"yellow\"), file.type=\"jpg\")\n";
        print GR "mlmMVP <- MVP(phe=pheno,geno=geno, map=map, K=Kinship, nPC.MLM=3, priority=\"speed\", ncpus=2, vc.method=\"BRENT\", maxLoop=10, method.bin=\"EMMA\", threshold=0.05, method=c(\"MLM\"))\n";
        print GR "MVP.Report(mlmMVP, plot.type=\"m\", LOG10=TRUE, ylim=NULL, threshold=c(1e-8,1e-6),threshold.lty=c(1,2), col=c(\"grey60\",\"grey30\"), threshold.lwd=c(1,1), threshold.col=c(\"black\",\"grey\"), amplify=TRUE, chr.den.col=c(\"darkgreen\", \"yellow\", \"red\"),bin.size=1e6,signal.col=c(\"red\",\"green\"), signal.cex=c(1,1),signal.pch=c(19,19),file.type=\"jpg\",memo=\"\",dpi=300)\n";

        print GR "glmMVP <- MVP(phe=pheno,geno=geno, map=map, K=Kinship, nPC.GLM=3, priority=\"speed\", ncpus=2, vc.method=\"BRENT\", maxLoop=10, method.bin=\"EMMA\", threshold=0.05, method=c(\"GLM\"))\n";
        print GR "MVP.Report(glmMVP, plot.type=\"m\", LOG10=TRUE, ylim=NULL, threshold=c(1e-8,1e-6),threshold.lty=c(1,2), col=c(\"grey60\",\"grey30\"), threshold.lwd=c(1,1), threshold.col=c(\"black\",\"grey\"), amplify=TRUE, chr.den.col=c(\"darkgreen\", \"yellow\", \"red\"),bin.size=1e6,signal.col=c(\"red\",\"green\"), signal.cex=c(1,1),signal.pch=c(19,19),file.type=\"jpg\",memo=\"\",dpi=300)\n";

        print GR "cpuMVP <- MVP(phe=pheno,geno=geno, map=map, K=Kinship, nPC.FarmCPU=3, priority=\"speed\", ncpus=2, vc.method=\"BRENT\", maxLoop=10, method.bin=\"EMMA\", threshold=0.05, method=c(\"FarmCPU\"))\n";
        print GR "MVP.Report(cpuMVP, plot.type=\"m\", LOG10=TRUE, ylim=NULL, threshold=c(1e-8,1e-6),threshold.lty=c(1,2), col=c(\"grey60\",\"grey30\"), threshold.lwd=c(1,1), threshold.col=c(\"black\",\"grey\"), amplify=TRUE, chr.den.col=c(\"darkgreen\", \"yellow\", \"red\"),bin.size=1e6,signal.col=c(\"red\",\"green\"), signal.cex=c(1,1),signal.pch=c(19,19),file.type=\"jpg\",memo=\"\",dpi=300)\n";
        close GR;
        print GH "less $pheno|awk '{print \$1\"\\t\"\$3}' >$mvp/MVP.txt\n";
        print GH "Rscript $mvp/gwas.R\n";
        close GH;
}

sub structure{
        my $structuredir="$odir/06.structure";mkdir $structuredir unless -d $structuredir;
        my $admixture="$structuredir/admixture"; mkdir $admixture    unless -d $admixture;
        my $faststructure="$structuredir/faststructure";mkdir $faststructure unless -d $faststructure;

        open  STRUCTURE1,">$admixture/structure_step1.sh" or die $!;
        print STRUCTURE1 "$Bin/bin/plink --noweb --file $datadir/result  --chr-set $setNum  --geno 0.05 --make-bed --out $admixture/structure.QC\n";
        close STRUCTURE1;
        open POP,"$pop" or die $!;
        my %ind_pop;my %pop_list;
        my $popnum=0;
        while(my $line=<POP>){
                chomp $line;
                my ($ind,$pop)=split/\t/,$line;
                $ind_pop{$ind}=$pop;
                if(! exists $pop_list{$pop}){
                        $pop_list{$pop}=$popnum;
                        $popnum++;
                        }
        }
        close POP;
        my $indnum=0;
        my $ind_pop_list="$admixture/structure.ip_list.txt";
        open IND_POP,">$ind_pop_list" or die $!;
        foreach my $ind(sort { $ind_pop{$a} cmp $ind_pop{$b} }keys %ind_pop){
                my $pop = $ind_pop{$ind};
                print IND_POP join("\t",$ind,$indnum,$pop,$pop_list{$pop}),"\n";
                $indnum++;
                }
        close IND_POP;
        for(my $i=2;$i<=6;$i++){
                open  STRUCTURE2,">$admixture/structure_$i.sh" or die $!;
                print STRUCTURE2 "$Bin/bin/admixture  -j8  --cv $admixture/structure.QC.bed $i|tee log${i}.out\n";
                close  STRUCTURE2;
        }
        open  STRUCTURE3,">$admixture/structure_step3.sh" or die $!;
        print STRUCTURE3 "This script must run locally\n";
        print STRUCTURE3 "num=`sh $Bin/bin/chooseK.sh $admixture`\n";
        print STRUCTURE3 "$Bin/bin/perl $Bin/bin/convert_structure_output_FORMAT.pl $ind_pop_list $admixture/structure.QC.\$num.Q $admixture\n";
        print STRUCTURE3 "$Bin/bin/distructLinux1.1 -M $popnum -N $indnum -d $Bin/bin/drawparams -p sorted.popq -i sorted.indivq -o sorted.ps -K \$num\n";##absolute path is not allowed
        print STRUCTURE3 "$Bin/bin/ps2pdf $admixture/sorted.ps $admixture/sorted.pdf\n";
        close STRUCTURE3
;

        open  STRUCTURE4,">$faststructure/faststructure.sh" or die $!;
        print STRUCTURE4 "$Bin/bin/plink  --file $datadir/result --chr-set $setNum --noweb --make-bed --out $faststructure/structure\n";

        for my $i(2..6){
                print STRUCTURE4 "$Bin/bin/python $Bin/bin/faststructure/structure.py -K $i --input $faststructure/structure --output $faststructure/structure --full --seed=100\n";
        }
        print STRUCTURE4 "$Bin/bin/python $Bin/bin/faststructure/chooseK.py --input=$faststructure/structure > $faststructure/chooseK.out\n";
        print STRUCTURE4 "best_k=`tail $faststructure/chooseK.out -n +2|cut -d  \" \" -f 10`\n";
        print STRUCTURE4 "$Bin/bin/python $Bin/bin/faststructure/distruct.py -K \$best_k --input=$faststructure/structure --title='' --output=$faststructure/fast_structure.svg\n";
        print STRUCTURE4 "convert -density 300 $faststructure/fast_structure.svg $faststructure/fast_structure.png\n";
        close STRUCTURE4;
}
