# This program does maximal dependence decomposition using a data file containing one nucleotide sequence each line
# output "answer.txt". 
# This program gives the position that weights most in the signal sequences. And the calc the new weight matrix that gives the weights of the remaining positions conditioned on this "max" position.

# load data
dseq1<-as.matrix(read.table("input.txt"))

# set up parameters for data
(length<-length(dseq1[1,]))
(number<-length(dseq1[,1]))

# make an array to store the number of transitions of each type (e.g. a->a, a->g, etc) between each pair of nucleotide positions in the signal. 
counts=array(dim=c(length,length,4,4))
counts[,,,]<-0

# Order the nucleotides as a,g,c,t.
# fill in the array "counts" using 3 nested loops to go through all the positions in all seqs
for(i in 1:number)
{
	for (j in 1:(length-1))
	{	# record counts for half of all possible transitions (fill in upper half of the array "counts")
		for (k in (j+1):length)
		{
			if (dseq1[i,j]=='a' & dseq1[i,k]=='a') 
			{	counts[j,k,1,1]<-counts[j,k,1,1]+1}
			if (dseq1[i,j]=='a' & dseq1[i,k]=='g') 
			{	counts[j,k,1,2]<-counts[j,k,1,2]+1}
			if (dseq1[i,j]=='a' & dseq1[i,k]=='c') 
			{	counts[j,k,1,3]<-counts[j,k,1,3]+1}
			if (dseq1[i,j]=='a' & dseq1[i,k]=='t') 
			{	counts[j,k,1,4]<-counts[j,k,1,4]+1}
			if (dseq1[i,j]=='g' & dseq1[i,k]=='a') 
			{	counts[j,k,2,1]<-counts[j,k,2,1]+1}
			if (dseq1[i,j]=='g' & dseq1[i,k]=='g') 
			{	counts[j,k,2,2]<-counts[j,k,2,2]+1}
			if (dseq1[i,j]=='g' & dseq1[i,k]=='c') 
			{	counts[j,k,2,3]<-counts[j,k,2,3]+1}
			if (dseq1[i,j]=='g' & dseq1[i,k]=='t') 
			{	counts[j,k,2,4]<-counts[j,k,2,4]+1}
			if (dseq1[i,j]=='c' & dseq1[i,k]=='a') 
			{	counts[j,k,3,1]<-counts[j,k,3,1]+1}
			if (dseq1[i,j]=='c' & dseq1[i,k]=='g') 
			{	counts[j,k,3,2]<-counts[j,k,3,2]+1}
			if (dseq1[i,j]=='c' & dseq1[i,k]=='c') 
			{	counts[j,k,3,3]<-counts[j,k,3,3]+1}
			if (dseq1[i,j]=='c' & dseq1[i,k]=='t') 
			{	counts[j,k,3,4]<-counts[j,k,3,4]+1}
			if (dseq1[i,j]=='t' & dseq1[i,k]=='a') 
			{	counts[j,k,4,1]<-counts[j,k,4,1]+1}
			if (dseq1[i,j]=='t' & dseq1[i,k]=='g') 
			{	counts[j,k,4,2]<-counts[j,k,4,2]+1}
			if (dseq1[i,j]=='t' & dseq1[i,k]=='c') 
			{	counts[j,k,4,3]<-counts[j,k,4,3]+1}
			if (dseq1[i,j]=='t' & dseq1[i,k]=='t') 
			{	counts[j,k,4,4]<-counts[j,k,4,4]+1}		
		}
	}
}

# fill in the lower half of "array" by mirroring
for (i in length:2)
{	for (j in (i-1):1)
	{	for (m in 1:4)
		{	for (n in 1:4)
			{	counts[i,j,m,n]<-counts[j,i,m,n]}
		}
	}
}


# calculate Chi-square values and store in a 10x10 array 
assocs<-array(dim=c(length,length))
assocs[,]<-0
# loop through a, g, c, t to calc Chi-sq for each pair of positions
# can also use the command chisq.test(counts[y1, y2,  ,  ])$statistic
for (i in 1:length)
{	for (j in 1:length)
	{	for (m in 1:4)
		{	for (n in 1:4)
			{	e<-sum(counts[i,j,m,])*sum(counts[i,j,,n])/sum(counts[i,j,,])
				chi<-((counts[i,j,m,n]-e)^2)/e
				assocs[i,j]<-assocs[i,j]+chi
			}
		}
	}
}


# set diagnals to 0
for (i in 1:length)
{	assocs[i,i]<-0}

# calc P-val matrix
pvals<-array(dim=c(length,length))
pvals[,]<-0
for (i in 1:(length-1))
{	for (j in (i+1):length)
	{	pvals[i,j]<-chisq.test(counts[i, j,  ,  ])$p.value}
}

# Find the nucleotide position i that affects the others most strongly by summing the rows of assocs. Call this position wmax. 
# Note we cannot deal with cases where two positions have equal weight
rmax<-sum(assocs[1,])
wmax<-1
for (i in 2:length)
{	if (sum(assocs[i,])>rmax)
	{	rmax<-sum(assocs[i,])
		wmax<-i
	}
}


# use an array "count" to record the number of occurrance of a given nt at a given position provided that position wmax is fixed with a,g,c,or t for each of the 4 tables.
count<-array(dim=c(4,length,4))
count[,,]<-0
# fill in the array "count" using 2 nested loops to go through all the positions in all seqs
for(i in 1:number)
{
	for (j in 1:length)
	{	# record counts for each nt at each position except wmax
		if (dseq1[i,wmax]=='a') # position wmax is 'a'
		{
			if (dseq1[i,j]=='a') 
			{	count[1,j,1]<-count[1,j,1]+1}
			if (dseq1[i,j]=='g') 
			{	count[2,j,1]<-count[2,j,1]+1}
			if (dseq1[i,j]=='c') 
			{	count[3,j,1]<-count[3,j,1]+1}
			if (dseq1[i,j]=='t') 
			{	count[4,j,1]<-count[4,j,1]+1}		
		}
		if (dseq1[i,wmax]=='g') # position wmax is 'g'
		{
			if (dseq1[i,j]=='a') 
			{	count[1,j,2]<-count[1,j,2]+1}
			if (dseq1[i,j]=='g') 
			{	count[2,j,2]<-count[2,j,2]+1}
			if (dseq1[i,j]=='c') 
			{	count[3,j,2]<-count[3,j,2]+1}
			if (dseq1[i,j]=='t') 
			{	count[4,j,2]<-count[4,j,2]+1}		
		}
		if (dseq1[i,wmax]=='c') # position wmax is 'c'
		{
			if (dseq1[i,j]=='a') 
			{	count[1,j,3]<-count[1,j,3]+1}
			if (dseq1[i,j]=='g') 
			{	count[2,j,3]<-count[2,j,3]+1}
			if (dseq1[i,j]=='c') 
			{	count[3,j,3]<-count[3,j,3]+1}
			if (dseq1[i,j]=='t') 
			{	count[4,j,3]<-count[4,j,3]+1}		
		}
		if (dseq1[i,wmax]=='t') # position wmax is 't'
		{
			if (dseq1[i,j]=='a') 
			{	count[1,j,4]<-count[1,j,4]+1}
			if (dseq1[i,j]=='g') 
			{	count[2,j,4]<-count[2,j,4]+1}
			if (dseq1[i,j]=='c') 
			{	count[3,j,4]<-count[3,j,4]+1}
			if (dseq1[i,j]=='t') 
			{	count[4,j,4]<-count[4,j,4]+1}		
		}
	}
}

# record px, the frequency of given nt at position wmax
px<-c(0,0,0,0)
px[1]<-sum(dseq1[,wmax]=="a")
px[2]<-sum(dseq1[,wmax]=="g")
px[3]<-sum(dseq1[,wmax]=="c")
px[4]<-sum(dseq1[,wmax]=="t")
(px<-px/number)

#make weight matrix by the condition j/wmax
weight1=matrix(nrow=4,ncol=length);
weight1[]<-0;
numA<-sum(dseq1[,wmax]=="a");
for (i in 1:number)
    {if (dseq1[i,wmax]=="a")
        {for (j in 1:length)
            {if (dseq1[i,j]=="a"){weight1[1,j]=weight1[1,j]+(1/numA)}
             if (dseq1[i,j]=="g"){weight1[2,j]=weight1[2,j]+(1/numA)}
             if (dseq1[i,j]=="c"){weight1[3,j]=weight1[3,j]+(1/numA)}
             if (dseq1[i,j]=="t"){weight1[4,j]=weight1[4,j]+(1/numA)}
            }
        }
    }

weight2=matrix(nrow=4,ncol=length);
weight2[]<-0;
numG<-sum(dseq1[,wmax]=="g");
for (i in 1:number)
    {if (dseq1[i,wmax]=="g")
        {for (j in 1:length)
            {if (dseq1[i,j]=="a"){weight2[1,j]=weight2[1,j]+(1/numG)}
             if (dseq1[i,j]=="g"){weight2[2,j]=weight2[2,j]+(1/numG)}
             if (dseq1[i,j]=="c"){weight2[3,j]=weight2[3,j]+(1/numG)}
             if (dseq1[i,j]=="t"){weight2[4,j]=weight2[4,j]+(1/numG)}
            }
        }
    }

weight3=matrix(nrow=4,ncol=length);
weight3[]<-0;
numC<-sum(dseq1[,wmax]=="c");
for (i in 1:number)
    {if (dseq1[i,wmax]=="c")
        {for (j in 1:length)
            {if (dseq1[i,j]=="a"){weight3[1,j]=weight3[1,j]+(1/numC)}
             if (dseq1[i,j]=="g"){weight3[2,j]=weight3[2,j]+(1/numC)}
             if (dseq1[i,j]=="c"){weight3[3,j]=weight3[3,j]+(1/numC)}
             if (dseq1[i,j]=="t"){weight3[4,j]=weight3[4,j]+(1/numC)}
            }
        }
    }

weight4=matrix(nrow=4,ncol=length);
weight4[]<-0;
numT<-sum(dseq1[,wmax]=="t");
for (i in 1:number)
    {if (dseq1[i,wmax]=="t")
        {for (j in 1:length)
            {if (dseq1[i,j]=="a"){weight4[1,j]=weight4[1,j]+(1/numT)}
             if (dseq1[i,j]=="g"){weight4[2,j]=weight4[2,j]+(1/numT)}
             if (dseq1[i,j]=="c"){weight4[3,j]=weight4[3,j]+(1/numT)}
             if (dseq1[i,j]=="t"){weight4[4,j]=weight4[4,j]+(1/numT)}
            }
        }
    }


# Given the model, calculate the probability of the signal being a specified sequence eg. seq1='gccattgtaa'.
seq1<-array(c('g','c','c','a','t','t','g','t','a','a'),dim=c(1,length))
pseq1<-1
if (seq1[wmax]=='a')
{
	for (i in 1:length)
	{	if (seq1[i]=='a' & i != wmax)
		{	pseq1<-pseq1*weight1[1,i]}
		if (seq1[i]=='g' & i != wmax)
		{	pseq1<-pseq1*weight1[2,i]}
		if (seq1[i]=='c' & i != wmax)
		{	pseq1<-pseq1*weight1[3,i]}
		if (seq1[i]=='t' & i != wmax)
		{	pseq1<-pseq1*weight1[4,i]}
	}
	pseq1
	pseq1<-pseq1*px[1]
}
if (seq1[wmax]=='g')
{
	for (i in 1:length)
	{	if (seq1[i]=='a' & i != wmax)
		{	pseq1<-pseq1*weight2[1,i]}
		if (seq1[i]=='g' & i != wmax)
		{	pseq1<-pseq1*weight2[2,i]}
		if (seq1[i]=='c' & i != wmax)
		{	pseq1<-pseq1*weight2[3,i]}
		if (seq1[i]=='t' & i != wmax)
		{	pseq1<-pseq1*weight2[4,i]}
	}
	pseq1
	pseq1<-pseq1*px[2]
}
if (seq1[wmax]=='c')
{
	for (i in 1:length)
	{	if (seq1[i]=='a' & i != wmax)
		{	pseq1<-pseq1*weight3[1,i]}
		if (seq1[i]=='g' & i != wmax)
		{	pseq1<-pseq1*weight3[2,i]}
		if (seq1[i]=='c' & i != wmax)
		{	pseq1<-pseq1*weight3[3,i]}
		if (seq1[i]=='t' & i != wmax)
		{	pseq1<-pseq1*weight3[4,i]}
	}
	pseq1
	pseq1<-pseq1*px[3]
}
if (seq1[wmax]=='t')
{
	for (i in 1:length)
	{	if (seq1[i]=='a' & i != wmax)
		{	pseq1<-pseq1*weight4[1,i]}
		if (seq1[i]=='g' & i != wmax)
		{	pseq1<-pseq1*weight4[2,i]}
		if (seq1[i]=='c' & i != wmax)
		{	pseq1<-pseq1*weight4[3,i]}
		if (seq1[i]=='t' & i != wmax)
		{	pseq1<-pseq1*weight4[4,i]}
	}
	pseq1
	pseq1<-pseq1*px[4]
}
	

# output results
# counts
print("the chi-squared values are:")
assocs
pvals
print("wmax=");wmax
print("weight1:")
weight1
print("weight2:")
weight2
print("weight3:")
weight3
print("weight4:")
weight4
print("the probability that the word would be sequence gccattgtaa is:")
pseq1


