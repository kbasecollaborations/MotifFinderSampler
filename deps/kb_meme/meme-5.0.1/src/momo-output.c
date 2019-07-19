#include "motif-in.h"
#include "config.h"
#include "momo-output.h"
#include "momo-html-string.h"
#include "io.h"
#include "motif-spec.h"
#include "seq-reader-from-fasta.h"
#include "ceqlogo.h"
#include "momo-algorithm.h"
#include "momo-modl.h"
#include "fisher_exact.h"

const int MAX_HTML_MATCHES = 1000;

/**********************************************************************
 * This function saves MOMO results as a MEME motif file.
 *********************************************************************/
void print_momo_text_file(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  fprintf(momo_file, "MEME version %s\n\n", VERSION);
  fprintf(momo_file, "Alphabet= %.*s\n\n", 20, summary.alph->symbols + 1);
  fprintf(momo_file, "Background letter frequencies\n");
  
  const char* alph_letters = summary.alph_letters;

  ARRAY_T * bg_freqs = summary.bg_freqs;
  
  int i;
  int j;
  int k;
  int l;
  
  for (i = 0; i < strlen(alph_letters); ++i) {
    fprintf(momo_file, "%c %8.6f ", alph_letters[i], get_array_item_defcheck(i, bg_freqs));
  }
  fprintf(momo_file, "\n\n");

  ARRAYLST_T * mod_table_keys = summary.mod_table_keys;
  
  for (i = 0; i < arraylst_size(mod_table_keys); i++) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    ARRAYLST_T* motifs = mod_entry->motifinfos;
    for (j = 0; j < arraylst_size(motifs); ++j) {
      MOTIF_INFO_T* motifinfo = arraylst_get(j, motifs);
      MOTIF_T* motif = motifinfo->motif;

      fprintf(momo_file, "MOTIF %s\n", get_motif_id(motif));
      unsigned long motif_list_size = (unsigned long) arraylst_size(motifinfo->seqs);
      if (motif->log_evalue > -BIG) {
	double m1, e1;
	exp10_logx(motif->log_evalue/log(10.0), m1, e1, 1);
	fprintf(momo_file, "letter-probability matrix: alength= %d w= %d nsites= %lu E= %3.1fe%+04.0f\n", (int) (strlen(alph_letters)), options.width, motif_list_size, m1, e1);
      } else {
	fprintf(momo_file, "letter-probability matrix: alength= %d w= %d nsites= %lu E= 0\n", (int) (strlen(alph_letters)), options.width, motif_list_size);
      }

      MATRIX_T* freqs = get_motif_freqs(motif);
      for (k = 0; k < options.width; k++) {
	for (l = 0; l < strlen(alph_letters); l++) {
	  fprintf(momo_file, "%8.6f\t", get_matrix_cell_defcheck(k, l, freqs));
	}
	fprintf(momo_file, "\n");
      }
      fprintf(momo_file, "\n");
    }
  }
};

void momo_print_version(FILE *momo_file) {
  fprintf(momo_file, "MoMo version %s (Release date: %s)", VERSION, ARCHIVE_DATE);
};

void momo_print_command_line(FILE *momo_file, MOMO_OPTIONS_T options) {
  fputs(options.command_line, momo_file);
};

void print_algorithm_name(MOMO_OPTIONS_T options, FILE *fp) {
  ALGORITHM_T algorithm = options.algorithm;
  if (algorithm == Simple) {
    fprintf(fp, "simple");
  } else if (algorithm == Motifx) {
    fprintf(fp, "motif-x");
  } else if (algorithm == Modl) { 
    fprintf(fp, "MoDL");
  } else {
    fprintf(fp, "UNKNOWN");
  }
}

void momo_print_parameters(FILE *momo_file, MOMO_OPTIONS_T options) {
  int i;
  
  fprintf(momo_file, "PARAMETERS:\n\n");
  
  // Algorithm
  fprintf(momo_file, "algorithm: "); print_algorithm_name(options, momo_file); fprintf(momo_file, "\n");
  
  // PTM files
  ARRAYLST_T* phospho_filenames = options.phospho_filenames;
  fprintf(momo_file, "post-translationally modified peptide filenames: \n");
  for (i = 0; i < arraylst_size(phospho_filenames); ++i) {
    char* phospho_filename = arraylst_get(i, phospho_filenames);
    fprintf(momo_file, "\tfile %d: %s\n", i+1, phospho_filename);
  }

  // PTM types and mod pep column.
  FILETYPE_T filetype = options.filetype;
  FILETYPE_T bg_filetype = options.bg_filetype;
  if (filetype == Psm) {
    fprintf(momo_file, "PTM filetype: %s\n", options.psm_type ? options.psm_type : "PSM");
    fprintf(momo_file, "modified peptide column: '%s'\n", options.sequence_column);
  } else if (filetype == Prealigned) {
    fprintf(momo_file, "PTM filetype: Raw\n");
  } else { // fasta
    fprintf(momo_file, "PTM filetype: FASTA\n");
  }

  // Context sequence
  if (options.protein_database_filename) fprintf(momo_file, "protein database filename: %s\n", options.protein_database_filename);
  if (options.protein_database_filename) {
    if (bg_filetype == Psm) {
      fprintf(momo_file, "protein database format: PSM (error!)\n");
    } else if (bg_filetype == Prealigned) {
      fprintf(momo_file, "protein database format: Raw\n");
    } else { // fasta
      fprintf(momo_file, "protein database format: FASTA\n");
    }
  }
  
  fprintf(momo_file, "motif width: %d\n", options.width);

  fprintf(momo_file, "filter: %s\n", (options.filter ? "true" : "false"));
  if (options.filter) {
    fprintf(momo_file, "\tfilter field: '%s'\n", options.filter_field);
    if (options.filter_type == Le) {
      fprintf(momo_file, "\tfilter type: <=\n");
    } else if (options.filter_type == Lt) {
      fprintf(momo_file, "\tfilter type: <\n");
    } else if (options.filter_type == Eq) {
      fprintf(momo_file, "\tfilter type: =\n");
    } else if (options.filter_type == Gt) {
      fprintf(momo_file, "\tfilter type: >\n");
    } else if (options.filter_type == Ge) {
      fprintf(momo_file, "\tfilter type: >=\n");
    } else {
      fprintf(momo_file, "\tfilter type is unknown!\n");
    }
    fprintf(momo_file, "\tfilter threshold: %g\n", options.filter_threshold);
  }

  fprintf(momo_file, "remove unknowns: %s\n", (options.remove_unknowns ? "true" : "false"));
  fprintf(momo_file, "eliminate repeats: %s\n", options.eliminate_repeat_width > 0 ? "true" : "false");
  if (options.eliminate_repeat_width > 0)
    fprintf(momo_file, "\teliminate repeat width: %d\n", options.eliminate_repeat_width);
  fprintf(momo_file, "min occurrences: %d\n", options.min_occurrences);
  fprintf(momo_file, "single motif per mass: %s\n", (options.single_motif_per_mass ? "true" : "false"));
  fprintf(momo_file, "hash fasta: %s\n", (options.hash_fasta ? "true" : "false"));
  fprintf(momo_file, "\thash fasta width: %d\n", options.hash_fasta_width);

  //fprintf(momo_file, "allow clobber: %s\n", (options.allow_clobber ? "true" : "false"));
  //fprintf(momo_file, "html path: %s\n", options.html_path);
  //fprintf(momo_file, "output dirname: %s\n", options.output_dirname);
  //fprintf(momo_file, "text path: %s\n", options.text_path);
  //fprintf(momo_file, "html filename: %s\n", options.HTML_FILENAME);
  //fprintf(momo_file, "text filename: %s\n", options.TEXT_FILENAME);

  if (options.algorithm == Motifx) {
    //fprintf(momo_file, "count threshold: %d\n", options.count_threshold);
    fprintf(momo_file, "score threshold: %g\n", options.score_threshold);
    fprintf(momo_file, "p-value calculations: %s\n", options.harvard ? "inaccurate (emulate original motif-x)" : "accurate");
  }

  if (options.algorithm == Modl) {
    fprintf(momo_file, "max motifs: %d\n", options.max_motifs);
    fprintf(momo_file, "max iterations: %d\n", options.max_iterations);
    fprintf(momo_file, "max no decrease iterations: %d\n", options.max_no_decrease);
  }
  
  // ints

  fprintf(momo_file, "\n");

  
};

static int n_tests(
  char *motifid,
  int w,
  int alen
) {
  int i;
  int len = strlen(motifid);				// characters in motif name
  int n = w-1;						// non-key,non-X positions in motif
  // Get the number of non-key,non-X positions in the motif.
  for (i=0; i<len; i++) if (motifid[i] == 'X') n--;
  // Add one position for the final tests that failed if n < w-1.
  if (n < w-1) n++;
  int n_tests = alen * n * (w - (n+1)/2.0);
  //printf("motif %s len %d w %d n %d n_tests %d\n", motifid, len, w, n, n_tests);
  return(n_tests);
} // n_tests

void momo_print_summary(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  fprintf(momo_file, "  <li>\n");
  fprintf(momo_file, "    SUMMARY:\n");
  fprintf(momo_file, "    <ul>\n");
  fprintf(momo_file, "      <li>Algorithm: "); print_algorithm_name(options, momo_file); fprintf(momo_file, "</li>\n");
  fprintf(momo_file, "      <li>Number of Mods: %lu</li>\n", summary.num_mod);
  fprintf(momo_file, "      <li>Number of Mod Types: %lu</li>\n", summary.num_modtype);
  fprintf(momo_file, "      <li>Number of Mods Passing Filters: %lu</li>\n", summary.num_mod_passing);
  fprintf(momo_file, "      <li>Number of Mod Types Passing Filters: %lu</li>\n", summary.num_modtype_passing);
  fprintf(momo_file, "    </ul><br>\n");
  fprintf(momo_file, "  </li>\n");
  
  ARRAYLST_T * mod_table_keys = summary.mod_table_keys;
  
  int i, j, k;
  for (i = 0; i < arraylst_size(mod_table_keys); ++i) {
    HASH_TABLE_ENTRY * hash_entry = arraylst_get(i, mod_table_keys);
    MOD_INFO_T * mod_entry = (MOD_INFO_T *) hash_get_entry_value(hash_entry);
    ARRAYLST_T * motifs = mod_entry->motifinfos;
    for (j = 0; j < arraylst_size(motifs); ++j) {
      MOTIF_INFO_T* currmotifinfo = arraylst_get(j, motifs);
      MOTIF_T * currmotif = currmotifinfo->motif;
      char* motifid = currmotif->id + 1;
      char* motif_name = mm_malloc(strlen(motifid) + 5);
      motif_name[0] = '\0';
      strncat(motif_name, motifid, strlen(motifid));
      strncat(motif_name, ".png", 4);
      
      char* motif_file = mm_malloc(strlen(options.output_dirname) + strlen(motifid) + 2);
      motif_file[0] = '\0';
      strncat(motif_file, options.output_dirname, strlen(options.output_dirname));
      strncat(motif_file, "/", 1);
      strncat(motif_file, motifid, strlen(motifid));
      CL_create1(currmotif, false, false, "MoMo", motif_file, false, true);
      
      if (options.algorithm == Motifx || options.algorithm == Modl) {
        double fold_enrich = (double) currmotifinfo->fg_match*currmotifinfo->bg_size / (currmotifinfo->bg_match*currmotifinfo->fg_size);
        double ntests = (options.algorithm == Modl) ? currmotifinfo->n_tests : n_tests(motifid, currmotif->length, strlen(summary.alph_letters));
        double log_n_tests = log(ntests);
        double log_pvalue = getLogFETPvalue(currmotifinfo->fg_match, currmotifinfo->fg_size, currmotifinfo->bg_match, currmotifinfo->bg_size, false);
        double log_evalue = currmotif->log_evalue = log_pvalue + log_n_tests;
        double m1, e1;
	exp10_logx(log_evalue/log(10.0), m1, e1, 1);
        if (options.printp) { 
          double log_norm_pvalue = LOGEV(log_n_tests, log_pvalue);
          double m2, e2;
	  exp10_logx(log_norm_pvalue/log(10.0), m2, e2, 1);
          printf("ntests: %.1f norm_p-value: %4.2fe%+04.0f E-value: %4.2fe%+04.0f\n", ntests, m2, e2, m1, e1);
        }
        fprintf(momo_file, "  <li>final_pattern: %s score: %.2f foreground_matches: %d foreground_size: %d bg_matches: %d bg_size: %d fold_increase: %.2f E-value: %3.1fe%+04.0f<br>\n", motifid, currmotifinfo->score, currmotifinfo->fg_match, currmotifinfo->fg_size, currmotifinfo->bg_match, currmotifinfo->bg_size, fold_enrich, m1, e1);
      } else if (options.algorithm == Simple) {
        currmotif->log_evalue = 0;
        fprintf(momo_file, "  <li>final_pattern: %s foreground_matches: %d<br>\n", motifid, currmotifinfo->fg_match);
      }
      
      fprintf(momo_file, "      <img src=\"%s\" alt=\"sequence logo of motif\"><br>\n", motif_name);
      fprintf(momo_file, 
        "      <button onclick=\"change_display('occ_%d_%d')\">Show/Hide Motif Occurrences</button>\n"
        "      <div style='display:none' id='occ_%d_%d'>\n"
        "        <pre class='console'>\n",
        i, j, i, j);
      for (k = 0; k < arraylst_size(currmotifinfo->seqs); ++k) {
        char* curr_motifinfo_seq = arraylst_get(k, currmotifinfo->seqs);
        fprintf(momo_file, "%s\n", curr_motifinfo_seq);
      }
      fprintf(momo_file, 
        "      </pre>\n"
        "    </div>\n"
        "    <br><br>\n"
        "  </li>\n"
      );
      
      // cleanup
      myfree(motif_name);
      myfree(motif_file);
    }
    // Print out MoDL log
    if (options.algorithm == Modl) {
      fprintf(momo_file,
        "  <li>\n"
        "    <button onclick=\"change_display('log_%d')\"><b>Show/Hide MoDL Log:</b></button>\n"
        "    <div style='display:none' id='log_%d'>\n"
        "      <ul>\n"
        "        <li><b>MoDL Log</b></li>\n",
        i, i);
      MATRIX_T* bg_freqs = NULL;
      bg_freqs = get_count_matrix(bg_freqs, mod_entry->bg_seq_list, NULL, &options, &summary);
      for (j = 0; j < options.width; ++j) {
        for (k = 0; k < strlen(summary.alph_letters); ++k) {
          set_matrix_cell_defcheck(j,k,get_matrix_cell_defcheck(j,k,bg_freqs)/arraylst_size(mod_entry->bg_seq_list),bg_freqs);
        }
      }
      
      ARRAYLST_T* modl_ops = mod_entry->modl_ops;
      ARRAYLST_T* temp_list = arraylst_create();
      double minDL = INFINITY;
      for (j = 0; j < arraylst_size(modl_ops); ++j) {
        MODL_STEP_T* step = arraylst_get(j, modl_ops);
        fprintf(momo_file, "      <li>STEP: %d, DL: %g<br>\n", j, step->score);
        do_step(step, temp_list, bg_freqs, options.max_motifs, &options, &summary, mod_entry);
        if (step->score < minDL) minDL = step->score;
        for (k = 0; k < arraylst_size(temp_list); ++k) {
          fprintf(momo_file, "        %s<br>\n", regexmotif_to_string(arraylst_get(k, temp_list), mod_entry, &summary, &options));
        }
        fprintf(momo_file, "      </li>\n");
      }
      fprintf(momo_file, 
        "      </ul>\n"
        "    </div><br>\n"
        "    <b>Final DL: %g</b>\n" 
        "  <br><br></li>\n",
        minDL);
      free_matrix(bg_freqs);
    } // modl
  } // i
  
};


/**********************************************************************
 * This function saves MOMO results as an HTML file
 *********************************************************************/
void print_momo_html_file(FILE *momo_file, MOMO_OPTIONS_T options, SUMMARY_T summary) {
  const int MAX_TAG_SIZE = 1000;
  int html_string_size = strlen(momo_html_string);
  int i = 0;
  for (i = 0; i < html_string_size; ++i) {
    if (momo_html_string[i] != '@') {
      fputc(momo_html_string[i], momo_file);
    }
    else {
      char buffer[MAX_TAG_SIZE];
        ++i;
      int j = 0;
      while (momo_html_string[i] != '@' && j < (MAX_TAG_SIZE - 1)) {
        buffer[j] = momo_html_string[i];
        ++j;
        ++i;
      }
      if (momo_html_string[i] != '@') {
        die("MoMo tag buffer length exceeded\n");
      }
      buffer[j] = '\0';
      if (strcmp("version", buffer) == 0) {
        momo_print_version(momo_file);
      }
      else if (strcmp("command_line", buffer) == 0) {
        momo_print_command_line(momo_file, options);
      }
      else if (strcmp("summary", buffer) == 0) {
        momo_print_summary(momo_file, options, summary);
      }
      else if (strcmp("parameters", buffer) == 0) {
        momo_print_parameters(momo_file, options);
      }
    }
  }
}

void create_directory(MOMO_OPTIONS_T options) {
  
  const bool PRINT_WARNINGS = false;
  if (create_output_directory(
                              options.output_dirname,
                              options.allow_clobber,
                              PRINT_WARNINGS
                              )
      ) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", options.output_dirname);
  }
}

/**********************************************************************
 * This function saves the MOMO results as a set of files in a
 * directory:
 *
 *   html_filename will be the name of the HTML output
 *   text_filename will be the name of the plain text output
 *
 * allow_clobber will determine whether or not existing files will
 * be overwritten.
 *********************************************************************/
void print_momo_results(MOMO_OPTIONS_T options, SUMMARY_T summary) {
  
  // Create directory for motifs to have a location
  create_directory(options);
  
  // Print HTML (must come first so it can set the E-values)
  FILE *momo_file = fopen(options.html_path, "w");
  if (!momo_file) {
    die("Couldn't open file %s for output.\n", options.html_path);
  }
  print_momo_html_file(momo_file, options, summary);
  fclose(momo_file);

  // Print plain text.
  momo_file = fopen(options.text_path, "w");
  if (!momo_file) {
    die("Couldn't open file %s for output.\n", options.text_path);
  }
  print_momo_text_file(momo_file, options, summary);
  fclose(momo_file);
  
}
