#require(devtools)
#devtools::install_github("kgori/sigfit", build_opts = c("--no-resave-data", "--no-manual"))
require(sigfit)

plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/sigfit/sandbox/"
# going through tutorial https://htmlpreview.github.io/?https://github.com/kgori/sigfit/blob/master/doc/sigfit_vignette.html#installation

############ simulating mutations (comments from tutorial) #########
#"First of all, we need some mutational signatures to fit to our data. The line below loads the 30 mutational signatures published in COSMIC (v2)."
data("cosmic_signatures_v2")
#Let’s use these signatures to simulate some mutation data. This code will generate 20,000 mutations from a 4:3:2:1 mixture of signatures 1, 3, 7 and 11.
# %*% is matrix multiplication. so you're multipying the probabilities by the cosmic signatures
set.seed(1)
# this gets probabilities of each 3mer mutation type from 4 signatures, each signature contributing with the ratios .4 .3 .2 .1 (using matrix multiplcation)  -- this is tricky but cool
probs <- c(0.4, 0.3, 0.2, 0.1) %*% as.matrix(cosmic_signatures_v2[c(1, 3, 7, 11), ])
# then you use these probabilities to sample 20000 mutations once 
mutations <- matrix(rmultinom(1, 20000, probs), nrow = 1)

colnames(mutations) <- colnames(cosmic_signatures_v2) # reanme by triplet type
head(mutations)
# this gives raw counts
plot_spectrum(mutations, name = "Simulated counts") # 

############# fit cosmic signatures ###############
# The arguments warmup and iter specify the number of ‘warm-up’ (or ‘burn-in’) iterations and the total number of sampling iterations (including warm-up iterations), respectively. We recommend that the value of warmup be between one-third and half the value of iter
# fit_signatures
# 
#By default, rstan::sampling uses iter = 2000 and warmup = floor(iter/2); going below these values is generally not advised.
#In general, one should run as many MCMC iterations (iter) as one’s computer and patience allow, with time and memory being the major constraints.
# so try to up the iters
# iters = as many as possible, min 2k
# The chains argument can be used to parallelise the sampling process by using multiple Markov chains, the number of which should never exceed the number of processors (the rstan::sampling function uses chains = 4 by default). To run multiple chains in parallel, the number of processors first needs to be set using options(mc.cores = parallel::detectCores()).
# chains = parallelize (not more than num of processors)
# iters = as many as possible, min 2k
# warmup (leave as default)
# The fit_signatures function can be used to fit the COSMIC signatures to the simulated counts as follows.
# this is with counts; just for one sample so far (?)
mcmc_samples_fit <- fit_signatures(counts = mutations, 
                                   signatures = cosmic_signatures_v2,
                                   iter = 2000, 
                                   warmup = 1000, 
                                   chains = 1, 
                                   seed = 1756)

####### here's where you can run other models. 

# this just took a couple seconds
# For details about further arguments of fit_signatures, including alternative signature models, mutational opportunities and exposure priors, see the Other functionalities section at the end of this vignette, or type ?fit_signatures to read the documentation.
# eventually will want to use the other models, etc. but wait to get there
# go through results
# Once we have the result of the MCMC sampling in mcmc_samples_fit, we can retrieve the estimated exposures from it using the retrieve_pars function
# # This returns a named list with three matrices, one containing the mean exposures, and the others containing the values corresponding to the lower and upper limits of the highest posterior density (HPD) interval for each signature exposure in each sample (HPD intervals are the Bayesian alternative to classical confidence intervals). 
exposures <- retrieve_pars(mcmc_samples_fit, 
                           par = "exposures",hpd_prob = .95) # 
#This gives you the mean exposures, and lower and upper limits of posterior density (HPD) interval for each signature exposure (bayesian alt to confidence interval); default is 95% but can change using hpd_prob for intervals
names(exposures)
exposures$mean # mean expsoures. can see 
# In this case, because there is only one sample, the exposures are returned as a numeric vector. For cases with multiple samples, exposures$mean is a matrix with one row per sample and one column per signature.
# The entire posterior distribution of the signature exposures and other model parameters in the mcmc_samples_fit object can be further explored by means of the functions provided by the rstan package. In addition, ShinyStan can be easily used in R for visual exploration of the MCMC samples.

############ plotting ###############

# the plot_spectrum function allows visualisation of both mutational catalogues and mutational signatures

# hese plots can be produced even for catalogues and signatures with arbitrary mutation types (e.g. indel or rearrangement signatures; see the Other functionalities section below).

# try arbitrary
# okay so you can plot this easily with any categories
test = tibble("A>T"=79,"C>T"=100,"A>G"=20) # made this up 
plot_spectrum(test, name = "Simulated counts") # 

# The path to an output PDF file can be specified using the pdf_path argument
# When plotting to a PDF, each function will automatically define appropriate graphical parameters, such as plot sizes and margins.
# custom plot sizes can be specified using the pdf_width and pdf_height arguments.
?plot_spectrum # for details

# plot_exposures makes barplots of exposures across the samle set.
# give it object resulting from MCMC sampling (mcmc_samples argument) 
# or exposures argument
plot_exposures(mcmc_samples=mcmc_samples_fit)
plot_exposures(counts=mutations,exposures=exposures) # could plot with counts and exposures
# The bars in this plot are coloured blue if the estimated exposure value is ‘sufficiently non-zero’, 
# In practice, ‘sufficiently non-zero’ means that the lower end of the Bayesian HPD interval (see the previous section) is above a threshold value close to zero (by default 0.01, and adjustable via the thresh argument). In this example, sigfit has identified the 4 signatures used to construct the sample.

########### rerun fit sigs , restricting to just those identified signatures #########
#  which have been highlighted as ‘sufficiently non-zero’ in the plot above, in order to obtain more-accurate estimates of exposures
selected_signatures <- c(1, 3, 7, 11) # based on 'blue' parts of plot

mcmc_samples_fit_2 <- fit_signatures(mutations, 
                                     # restricting just to selected signatures
                                     cosmic_signatures_v2[selected_signatures, ],
                                     iter = 2000, 
                                     warmup = 1000, 
                                     chains = 1, 
                                     seed = 1756)
mcmc_samples_fit_2
exposures2_restricted <- retrieve_pars(mcmc_samples_fit_2,par="exposures")
exposures2_restricted # refine exposures

######### plot reconstruction #########
plot_reconstruction(mcmc_samples=mcmc_samples_fit)
# this will show contributions of all the possible signatures
# whereas this will show it with the ones I restricted to:
plot_reconstruction(mcmc_samples=mcmc_samples_fit_2,pdf_path = paste0(plotdir,"test1.pdf")) #
# I like these plots much more --- can help distuingish the CpG>T from other signautres I bet!
# can use plot_all to plot spectrum, exposures and reconstructions
# can only do this to a pdf file 
plot_all(mcmc_samples=mcmc_samples_fit,out_path = plotdir,prefix="test2")

########### Extracting signatures from multiple samples ############
data("variants_21breast")
head(variants_21breast)
# okay so this is a way to make a variant table (but I'm going to use the spectrum type later)
dim(variants_21breast)
# note: ariants which are found in more than one sample need to be included multiple times in the table, using different sample IDs. 
#  create catalogues:
counts_21breast <- build_catalogues(variants_21breast) # only works for triplets
# but I dont need build_catalogues() -- going to make my own catalogues.
head(counts_21breast)
# so to format my data for this, I will want
# samples as rows, mtuation types as column names (names of mutaiton types can be formatted any way I think)
# The mutational catalogues are stored as a matrix of mutation counts,
# , where each row refers to a sample and each column corresponds to a trinucleotide mutation type. (This particular set of 21 catalogues is also available via data("counts_21breast").)
# Note that the build_catalogues function only admits mutational catalogues defined over the typical 96 trinucleotide mutation types (or 192 types if using transcriptional strand information). Mutational catalogues defined over other sets of mutation types need to be generated by the user.

counts_21breast[1:5, 1:8]
# For tables containing multiple catalogues, this function will produce one plot per catalogue, 
# making the use of an output PDF file (pdf_path argument) more convenient.
#par(mfrow = c(7, 3)) # having trouble getting the grid to work
plot_spectrum(counts_21breast,pdf_path = paste0(plotdir,"multisample.pdf"))
# this plots them all in one file

###### TO DO: make sure I figure out 'opportunities' even if tutorial doesn't
######### YOU ARE HERE

############## de novo signature extraction ##############
# need to use nsignatures (single or multi)
# # he recommended approach is to first run the function for a small number of iterations and a reasonably wide range of numbers of signatures (e.g. nsignatures = 2:8). When nsignatures is a range of values, sigfit will automatically determine the most plausible number of signatures present in the data (which is done by assessing goodness-of-fit through the calculate_gof and plot_gof functions).
# this is slower 

# Note that, for the correct determination of the number of signatures, nsignatures should include at least four values, and should extend beyond the ‘reasonable’ range for the expected number of signatures (e.g. if the number of signatures is expected to be between 3 and 8, nsignatures could take values 2:10). Also, note that for ranges of nsignatures the output will be a list, in which element [[N]] corresponds to the extraction results for nsignatures = N.


# got an error when running with 2000 iterations; needed to up to 
# The number of Markov chains in extract_signatures has a default value of chains = 1. We do not recommend the use of multiple chains for signature extraction, as this can result in a problem known as ‘label switching’ which can invalidate the inferences.
# you can plot goodness of fit;
#in which element [[N]] corresponds to the extraction results for nsignatures = N.
mcmc_samples_extr <- extract_signatures(counts = counts_21breast,
                                        nsignatures = 2:7,
                                        iter = 1000, 
                                        seed = 1756)
# yields warnings that posterios means and medians may be unreliable 
# hm

######## plot goodness of fit ########
png(paste0(plotdir,"goodnessoffit.png"),res = 300,height = 5,width=7)
par(mar=c(1,1,1,1))
best = plot_gof(mcmc_samples_extr) # this will save 'best' as a number
dev.off()
# okay plot manually:
#ggplot(data.frame(nS=mcmc_samples_extr$gof$nS,gof=mcmc_samples_extr$gof$gof),aes(x=nS,y=gof))+
#  geom_point()+
#  geom_line()

  # okay so you can see in this plot the best number of signatures (elbow of plot)
# how to get one particular signature:
# okay note that 'best' is apparently the index?
#DONT DO THIS: best=mcmc_samples_extr$gof$best # gets number of 'best' signature count # okay no this is maybe the best 'index'? don't use this. because here is saying best is 3 when it's 4. so maybe it's position 3 in the plot ? but not the actual best
# this part of the program isnt' great


# Next, it is recommended to re-run the extract_signatures function, this time with nsignatures = 4 and a much larger number of iterations, to obtain more-accurate estimates. We will skip this step in the present example.

# so rerun 
######### THIS IS SLOW: don't rerun ############
mcmc_samples_extr_best_LOTSOFITER <- extract_signatures(counts = counts_21breast,
                                        nsignatures = best,
                                        iter = 10000, 
                                        seed = 1756)

########## retrieving and extracting signatures ###########
# retrieve_pars function with par = "signatures" or par = "exposures"
sigs <- retrieve_pars(mcmc_samples_extr[[best]],par="signatures")
exposures <- retrieve_pars(mcmc_samples_extr[[best]],par="exposures")

# you can see which Cosmic signatures are closest to these sigs:
match_signatures(sigs, cosmic_signatures_v2) # how to save this output?

# this gives an assignment
# In this case, the closest matches of the extracted signatures A, B, C and D are, respectively, COSMIC signatures 2, 3, 13 and 1.
# As a last remark, once that more-accurate signatures have been extracted by re-running extract_signatures for a single value of nsignatures (as recommended above), it might be useful to re-fit these signatures back to the original catalogues, as this sometimes results in exposure estimates that are more accurate than the ones obtained by signature extraction. 
#Assuming that the MCMC samples from the definitive signature extraction [ with lots of iter ] are stored in an object named mcmc_samples_extr_LOTSOFITER, re-fitting is done as follows.

signatures <- retrieve_pars(mcmc_samples_extr_LOTSOFITER, "signatures")
mcmc_samples_refit <- fit_signatures(counts = counts_21breast,
                                     signatures = signatures,
                                     iter = 2000,
                                     warmup = 1000)
exposures <- retrieve_pars(mcmc_samples_refit, "exposures") # refine exposures

########### alternate models ############
# fit_signatures and extract_signatures  can use
# Poisson, negbin, normal
# I want to try default (nmf); oisson and negbin

############ using mutational opportunities (targets) ###########
#  For catalogues composed of 96 trinucleotide mutation types, these opportunities are given by the frequency of each trinucleotide in the target sequence (usually a genome or exome). 
# Mutational opportunities can be specified using the opportunities argument in the functions fit_signatures and extract_signatures. Opportunities can be provided as a numeric matrix with one row per sample and one column per mutation type. 
# # Importantly, using mutational opportunities during signature extraction will result in signatures which are already normalised by the given opportunities (that is, their mutation probabilities are not relative to the composition of the target sequence). In contrast, signatures extracted without using opportunities are relative to the sequence composition; this is not a problem unless the signatures are to be subsequently applied to other types of sequences. (See below for details on how to convert signatures between different sets of opportunities.)
# priors for signatures and exposures via the arguments sig_prior and exp_prior. These priors allow the incorporation of a priori knowledge about these parameters; for instance, using certain signatures from COSMIC as our signature priors reflects our belief that the signatures found in the data should be very similar to these. Given their potential for introducing subjective biases into the analysis, the use of priors should be strictly limited to situations where there is both strong evidence and good reason to use them. Priors must be numeric matrices with the same dimensions as the signature and exposure matrices to be inferred: signature priors must have one row per signature and one column per mutation type, while exposure priors must have one row per sample and one column per signature.
#  if you use appropriate opportunities, you don't need to convert -- signatures will not be relative to oppos. ut if you don't then you may need to convert usign convert_signatures (used from genome > exome) but can also convert between two vectors of opportunities 
# between any two given vectors of opportunities (with one element per mutation type).
# be careful if converting !! lots to do. better to just use opportunities?

######### discover rare signatures : Fit -Ext ##############
# One novelty in sigfit is the introduction of Fit-Ext models, which are able to extract mutational signatures de novo while fitting a set of predefined signatures that are already known to be present in the data. Such models are useful for the discovery of rare signatures for which there is some qualitative evidence, but insufficient support as to enable deconvolution via conventional methods, or in cases where only signature fitting is possible, yet the data clearly display a mutational pattern which is not captured by any of the available signatures.
# this is what I want to do with aging signatures!
# fit_extract_signatures
# a matrix of known signatures to be fitted needs to be provided via the signatures argument (as in fit_signatures), and that the number of additional signatures to extract is specified using the num_extra_sigs argument
# can only be a signle integer (not a range)
?fit_extract_signatures
# trying it with having cosmic 1 be fixed 
test_fix_ext <- fit_extract_signatures(mutations,signatures = cosmic_signatures_v2[c(1, 3, 7, 11),],num_extra_sigs = 2) # trying to first fit the ones I knew were part of this, then extract two more
# how to assess fit 

fix_ext_sig <- retrieve_pars(test_fix_ext,"signatures")
# first set will be the ones we defined, plus two extra
fix_ext_exp <- retrieve_pars(test_fix_ext,"exposures")
match_signatures(fix_ext_sig,cosmic_signatures_v2) # okay so it finds two more signs but very low exposures 
plot_exposures(counts=mutations,exposures=fix_ext_exp) # could plot with counts and exposures
# shows that their exposure isn't significant -- so that is cool! 
plot_reconstruction(mcmc_samples=test_fix_ext)

##################### okay let's try it ####################
agingSignatures <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/sigfit/jonsson.agingsignatures.tables9.convertedtoproportionsinexcel.txt",header=T,sep="\t",row.names = "signature") # make a table of these with signature name as rowname (so it's like the cosmic format)


######## get my spectra ###############
spectradir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS/20211216.plots.projectedCounts.epsilon.1.sepCpGs.no.addingApes_maskedMouseVaq.addingMultinomialDownsampling.MorePlotsForPAG.newTimeTree.MantelTest/"

spectra=read.table(paste0(spectradir,"AllProcessedSpectra.AllConditions.AllVariables.AllSpecies.txt"),header=T)
# need to get counts and targets # 
# get from sript 

allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/20211123_All_SpectrumCounts_ProjectedDownTo.10.Haploids.includesApes.RepMaskedMice.txt",header=T)
# need to write out downsampled counts and targets
# just get one to start with? write out in sigfit fmt
# want downsampled counts -- are these in here? 
