function buildMotifobj(searchmotifs, basedict) {
	var returndict = {};
	const strands = Object.entries(basedict);
	for(const [feat_id, motifs] of strands) {
		const matches = Object.entries(motifs);
		for(const [motif, locations] of matches) {
			for(i=0;i<searchmotifs.length;i++) {
				if(motif == searchmotifs[i]) {
					if (typeof returndict[feat_id] == 'undefined') {
						returndict[feat_id] = {};
						returndict[feat_id][motif] = [];

						for(q=0;q<locations.length;q++) {
							returndict[feat_id][motif].push(locations[q]);
						}
					} else if (typeof returndict[feat_id][motif] == 'undefined') {
						returndict[feat_id][motif] = [];
						// each sequence entry is [start, end, orientation]
						for(q=0;q<locations.length;q++) {
							returndict[feat_id][motif].push(locations[q]);
						}
					} else {
						for(q=0;q<locations.length;q++) {
							returndict[feat_id][motif].push(locations[q]);
						}
					}
				}
			}
		}
	}
	return returndict;
}

var options = {
  "ymax": 2,
	"ylabel": "bits",
	"showaxis": false
 };

var loc_dict = {};
var first_motif = {};
var first_motif_seq;
var loc_len = 0;

$("#seq_loc_load").css('display','none');

fetch('ReportMotif.json')
	.then(function(response) {
		return response.json();
	})
	.then(function(motifset) {
		motifs = motifset["Motifs"];

		$(motifs).each(function(i) {
			var locs = this['Motif_Locations'];
			var canvaselem = "logo_canvas_from_json_"+(i+1);
			var revcompcanvaselem = "logo_canvas_from_json_revcomp_"+(i+1);
			var chkbox = (i==0) ? (' checked') : ('');

			if(i == 0) {
				first_motif_seq = this['Iupac_sequence'];
			}

			var inserthtml = 	'<div class="motif-grid-item"><input type="checkbox" name="chk_'+this['Iupac_sequence']+'" onchange="triggerChk();" value="'+this['Iupac_sequence']+'" '+chkbox+'></div>'+
												'<div class="motif-grid-item">'+(i+1)+'</div>'+
												'<div class="motif-grid-item">'+this['Iupac_sequence']+'</div>'+
												'<div class="motif-grid-logo"><canvas id="'+canvaselem+'"></canvas></div>'+
												'<div class="motif-grid-logo"><canvas id="'+revcompcanvaselem+'"></canvas></div>'+
												'<div class="motif-grid-item"><a class="download-char" href="#"><span>&#8628;</span></a></div>';

			$("#motif-grid").append(inserthtml);

			for(j=0;j<locs.length;j++) {
				feat_id = locs[j]['sequence_id'];
				seq = this['Iupac_sequence'];

				if (typeof loc_dict[feat_id] == 'undefined') {
					loc_dict[feat_id] = {};
					loc_dict[feat_id][seq] = [];
					loc_dict[feat_id][seq].push([locs[j]['start'], locs[j]['end'], locs[j]['orientation']]);
				} else if (typeof loc_dict[feat_id][seq] == 'undefined') {
					loc_dict[feat_id][seq] = [];
					// each sequence entry is [start, end, orientation]
					loc_dict[feat_id][seq].push([locs[j]['start'], locs[j]['end'], locs[j]['orientation']]);
				} else {
					loc_dict[feat_id][seq].push([locs[j]['start'], locs[j]['end'], locs[j]['orientation']]);
				}
			}

			sequence_logo(document.getElementById(canvaselem), 300, 100, PPM2PWM(this['PWM']), options);
			sequence_logo(document.getElementById(revcompcanvaselem), 300, 100, revcomplogomatrix(PPM2PWM(this['PWM'])), options);

			PPM2PWM(this['PWM']);
		});

		loc_len = Object.keys(loc_dict).length;
		seq_location(document.getElementById("seq_loc"), 1000, 100, buildMotifobj([first_motif_seq], loc_dict), options);
		$("input[name='chk_"+first_motif_seq+"']").parent().css('background', window.motifcolors[first_motif_seq]);
	});

	function triggerChk() {
		$("#seq_loc").css('display','none');
		$("#seq_loc_load").css('display', 'inline');
		var listofchkd = [];

		$("div.motif-grid-item input[type='checkbox']:checked").each(function () {
			listofchkd.push($(this).val());
		});

		seq_location(document.getElementById("seq_loc"), 1000, 100, buildMotifobj(listofchkd, loc_dict), options);

		$('div.motif-grid-item input[type="checkbox"]').each(function () {
			if($(this).is(":checked")) {
				$(this).parent().css('background', window.motifcolors[$(this).val()]);
			} else {
				$(this).parent().css('background', 'none');
			}
		});

		$("#seq_loc_load").css('display', 'none');
		$("#seq_loc").css('display','inline');
	}
