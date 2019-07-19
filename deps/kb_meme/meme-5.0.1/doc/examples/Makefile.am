##
## Process this file with automake to produce Makefile.in
##

##############################################

AUTOMAKE_OPTIONS = no-dependencies


SUBDIRS = \
  compute_prior_dist_example_output_files \
  sample_opal_scripts

EXAMPLES = \
  sample-dna-aln.glam2 \
  sample-protein-aln.glam2 \
  sample-dna-motif.meme \
  sample-protein-motif.meme \
  small-tide-mod.txt

# load the EXAMPLE_OUTPUT_FILES variable
include examples.mk

examplesdir = ${memedocdir}/examples
export examplesdir
nobase_dist_examples_DATA = $(EXAMPLES) $(EXAMPLE_OUTPUT_FILES)

EXTRA_DIST = examples.mk $(EXAMPLES) $(EXAMPLE_OUTPUT_FILES)
