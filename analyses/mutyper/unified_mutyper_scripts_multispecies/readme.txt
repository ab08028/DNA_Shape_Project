########## readme for mutyper variants script ###########
# the issue: each dataset needs slightly different settings
# but i want to keep track of what files go in and be able to change them (eg masks, ancestral fastas) easily
# and have that be consistent for downstream scripts like mutyper targets as well
# so I've made config files for each species that load in their particular options
# should be smoother once processed by mutyper variants # downstream scripts may not need as many special options

# examples of special options
# some species need bad qual inds or relatives or subspecies removed (fw, gorilla, chimp,vaquita)
# some species need PASS sites selected, others need some nonPASS sites left in (fw)
# humans use --strict option but no one else does
# to keep track of everything can look at config files
