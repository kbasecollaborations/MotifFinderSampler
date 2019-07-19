
var alphabet = null;
var sequences = null;
var control_sequences = null;
var motifs = null;
var background = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "alphabet") {
    alphabet = controler;
    element.addEventListener("alphabet_changed", function (e) {
      if (sequences != null) sequences.set_custom_alphabet(e.detail.has_custom, e.detail.alphabet);
      if (control_sequences != null) control_sequences.set_custom_alphabet(e.detail.has_custom, e.detail.alphabet);
    }, false);
  } else if (id == "sequences") {
    sequences = controler;
    if (alphabet != null) {
      sequences.set_custom_alphabet(alphabet.has_custom_alphabet(), alphabet.get_custom_alphabet());
    }
  } else if (id == "control_sequences") {
    control_sequences = controler;
    if (alphabet != null) {
      control_sequences.set_custom_alphabet(alphabet.has_custom_alphabet(), alphabet.get_custom_alphabet());
    }
  } else if (id == "motifs") {
    motifs = controler;
  } else if (id == "background") {
    background = controler;
  }
}

function check() {
  "use strict";
  var alphs = null;
  if (alphabet != null) {
    alphs = alphabet.get_alphabets();
  }
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    if (alphs == null) alphs = sequences.get_alphabets();
  }
  if ($("generate_off").checked) {
    if (control_sequences != null) {
      if(!control_sequences.check(alphs)) return false;
      if (alphs == null) alphs = control_sequences.get_alphabets();
    }
  }
  if (motifs != null) {
    if(!motifs.check(alphs, alphabet != null && alphabet.get_alphabets() != null)) return false;
    if (alphs == null) alphs = motifs.get_alphabets();
  }
  if (!check_num_value("motif match log-odds fraction", "hit_lo_fraction", 0, 1, 0.25)) return false;
  if (!check_num_value("E-value report threshold", "evalue_report_threshold", 0, 1e300, 10)) return false;
  if (!check_job_details()) return false;
  if (background != null && !background.check(alphs)) return false;
  return true;
}

function options_changed() {
  if ($("scoring").value != "avg") return true;
  if ($("hit_lo_fraction").value != 0.25) return true;
  if ($("method").value != "fisher") return true;
  if ($("kmer").value != 2) return true;
  if (!/^\s*10\s*$/.test($("evalue_report_threshold").value)) return true;
  if (background != null && background.changed()) return true;
  return false;
}

function options_reset(evt) {
  //reset options
  $("scoring").value = "avg";
  $("hit_lo_fraction").value = 0.25;
  //$('hits_area').style.display = 'none';
  $("method").value = "fisher";
  $("kmer").value = 2;
  //$('fisher_pwm_area').style.display = 'none';
  $("evalue_report_threshold").value = 10;
  if (background != null) background.reset();
}

function fix_reset() {
  on_ch_generate();
  on_ch_scoring();
}

function on_form_submit(evt) {
  if (!check()) {
    evt.preventDefault();
  }
}

function on_form_reset(evt) {
  window.setTimeout(function(evt) {
    fix_reset();
  }, 50);
}

function on_ch_generate() {
  var method_value = $('method').value;
  $('control').style.display = $('generate_off').checked ? 'block' : 'none';
  // Only allow setting fisher pwm threshold if there is no control sequence file.
  //$('fisher_pwm_area').style.display = $('generate_none').checked && /(fisher)/.test(method_value) ? 'block' : 'none';
  // Only allow setting kmer value if using a shuffled control.
  $('kmer_area').style.display = $('generate_on').checked ? 'block' : 'none';
  // Allow a control set only with fisher and ranksum methods.
  if (! $('generate_none').checked && (! (/(fisher|ranksum)/.test(method_value)))) {
    $('method').selectIndex = 1;
    $('method').value = 'fisher';
  }
}

function on_ch_method() {
  var method_value = $('method').value;
  // Only allow 'totalhits' with 3dmhg and 4dmhg methods
  if (/(3dmhg|4dmhg)/.test(method_value)) {
    $('scoring').selectIndex = 4;
    $('scoring').value = 'totalhits';
  }
  // Allow a control set only with fisher and ranksum methods.
  if (! (/(fisher|ranksum)/.test(method_value))) {
    $('generate_none').checked = true;
  }
  // Only allow setting fisher pwm threshold if there is no control sequence file.
  //$('fisher_pwm_area').style.display = $('generate_none').checked && /(fisher)/.test(method_value) ? 'block' : 'none';
  // Only allow setting kmer value if using a shuffled control.
  $('kmer_area').style.display = $('generate_on').checked ? 'block' : 'none';
}

function on_ch_scoring() {
  var scoring_value = $('scoring').value;
  //$('hits_area').style.display = /(totalhits)/.test(scoring_value) ? 'block' : 'none';
  var method_value = $('method').value;
  if (! /(totalhits)/.test(scoring_value) && /(3dmhg|4dmhg)/.test(method_value)) {
    $('method').selectIndex = 1;
    $('method').value = 'fisher';
  }
}

function on_pageshow() {
  alphabet._radio_update(alphabet);
  motifs._source_update();
  sequences._source_update();
  control_sequences._source_update();
  background._source_update();
  update_compare();
}

function on_load() {
  // add listeners for the control generation mode
  $("generate_off").addEventListener("click", on_ch_generate, false);
  $("generate_on").addEventListener("click", on_ch_generate, false);
  $("generate_none").addEventListener("click", on_ch_generate, false);
  // add listener to the form to check the fields before submit
  $("ame_form").addEventListener("submit", on_form_submit, false);
  $("ame_form").addEventListener("reset", on_form_reset, false);
  $("method").addEventListener("change", on_ch_method, false);
  $("scoring").addEventListener("change", on_ch_scoring, false);
  window.addEventListener('pageshow', on_pageshow, false);
  on_ch_scoring();
  on_ch_method();
}

// add a load
(function() {
  "use strict";
  window.addEventListener("load", function load(evt) {
    "use strict";
    window.removeEventListener("load", load, false);
    on_load();
  }, false);
})();
