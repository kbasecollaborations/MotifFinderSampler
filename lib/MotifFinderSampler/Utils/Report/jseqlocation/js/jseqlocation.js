// from https://github.com/frogcat/canvas-arrow
(function(target) {
  if (!target || !target.prototype)
    return;
  target.prototype.arrow = function(startX, startY, endX, endY, controlPoints) {
    var dx = endX - startX;
    var dy = endY - startY;
    var len = Math.sqrt(dx * dx + dy * dy);
    var sin = dy / len;
    var cos = dx / len;
    var a = [];
    a.push(0, 0);
    for (var i = 0; i < controlPoints.length; i += 2) {
      var x = controlPoints[i];
      var y = controlPoints[i + 1];
      a.push(x < 0 ? len + x : x, y);
    }
    a.push(len, 0);
    for (var i = controlPoints.length; i > 0; i -= 2) {
      var x = controlPoints[i - 2];
      var y = controlPoints[i - 1];
      a.push(x < 0 ? len + x : x, -y);
    }
    a.push(0, 0);
    for (var i = 0; i < a.length; i += 2) {
      var x = a[i] * cos - a[i + 1] * sin + startX;
      var y = a[i] * sin + a[i + 1] * cos + startY;
      if (i === 0) this.moveTo(x, y);
      else this.lineTo(x, y);
    }
  };
})(CanvasRenderingContext2D);

function randomcolor() {
	return '#'+(Math.random()*0xFFFFFF<<0).toString(16);
}

function seq_location(element, width, height, seqs, options = null) {
	// Get setttings from options or default values.
	var defaults = {
		"colors": {"+": "green", "-": "blue", },
		"textcolor": "black",
		"bgcolor": "white",
		"padding": 6,
		"legendheight": 150,
		"lineheight": "16px",
		"labelfont": "12px \"Lucida Console\", Monaco, monospace",
		"legendfont": "14px \"Lucida Console\", Monaco, monospace",
		"linecolor": "black"
	};
	
	var settings = {};
	for(var p in defaults) {
		settings[p] = (options[p] == null) ? defaults[p] : options[p];
	}
	
	if(typeof window.motifcolors == 'undefined') {
		window.motifcolors = {};
	}
	
	// stat collection
	var upstreamstrt = 0;
	var downstreamend = 0;
	const strands = Object.entries(seqs);
	for(const [feat_id, motifs] of strands) {
		const matches = Object.entries(motifs);
		for(const [motif, locations] of matches) {
			for(i=0;i<locations.length;i++) {
				// 0: start, 1:end, 2:orientation
				if(locations[i][0] < locations[i][1]) { // if start < end 5'=> 3'					
					if(locations[i][0] < upstreamstrt) {
						upstreamstrt = locations[i][0];
					}
					if(locations[i][1] > downstreamend) {
						downstreamend = locations[i][1];
					}
				} else { // but if end < start 3'->5'
					if(locations[i][1] < upstreamstrt) {
						upstreamstrt = locations[i][1];
					}
					if(locations[i][0] > downstreamend) {
						downstreamend = locations[i][0];
					}
				}
			}
		}
	}
	
	var streamdelta = downstreamend - (upstreamstrt);
	upstreamstrt = upstreamstrt - Math.floor(0.05*streamdelta);
	downstreamend = downstreamend + Math.floor(0.05*streamdelta);
	streamdelta = downstreamend - (upstreamstrt);
	
	// box model calculations
	textboxheight = height - (2*settings.padding);
	textboxwidth = 0.2 * width - (2*settings.padding);	
	strandboxheight = height - (2*settings.padding);
	strandboxwidth = (0.8*width)-((2*settings.padding)+20);
	sendsestrandheight = Math.floor(strandboxheight/2);
	pixel2strand = strandboxwidth/streamdelta;
	
	// setup canvas element
	var legdfontsize = Number((settings.legendfont).substring(0,2));
	var legendheight = (Object.keys(window.motifcolors).length == 0) ? legdfontsize + (2*settings.padding) : (Math.ceil((Object.keys(window.motifcolors).length+1)/2) * legdfontsize) + (2*settings.padding) + (2*Math.floor(2%Object.keys(window.motifcolors).length));
	var legendwidth = width-(2*settings.padding);
	var legendstartx = 0+settings.padding;
	var legendstarty = (Object.keys(seqs).length) * height + settings.padding;
	totalheight = (Object.keys(seqs).length) * height + settings.legendheight;
	element.width = width;
	element.height = totalheight;
	var ctx = element.getContext("2d");
	ctx.save();
	seqnum = 0;
	
	for(const [feat_id, motifs] of strands) {
		// draw sequence location labels
		ctx.fillStyle="black";
		ctx.font = settings.labelfont;
		var p = settings.padding*(2*seqnum+1);
		var wraptextx = (p+(textboxheight*(seqnum+1)))-(textboxheight/2);
		ctx.fillText(feat_id, settings.padding, wraptextx);
		ctx.save();
				
		// draw strand diagram line and sense labels		
		ctx.font = "16px \"Lucida Console\", Monaco, monospace";		
		ctx.fillText("+",Math.floor(((0.2*width)+settings.padding)),Math.floor(((height/2)-8)+(height*seqnum)));
		ctx.fillText("-",Math.floor(((0.2*width)+settings.padding)),Math.floor(((height/2)+14)+(height*seqnum)));
		ctx.save();
		ctx.beginPath();
		ctx.moveTo(Math.floor(0.2*width),Math.floor((height/2)+(height*seqnum)));
		ctx.lineTo(width,Math.floor((height/2)+(height*seqnum)));
		ctx.lineWidth = 1;
		ctx.strokeStyle = settings.linecolor;
		ctx.stroke();
		
		const matches = Object.entries(motifs);
		for(const [motif, locations] of matches) {
			for(i=0;i<locations.length;i++) {
				/* 
					Calculations for sequence motif box models.
					pixel2strand - ratio of pixels to nucleotides on the current strand model -a
					px2start - pixels from beginning of strand diagram to start of motifbox
					pxofstrand - pixel length of the motif box
				*/
				var seqboxy = (locations[i][2] == "+") ? ((height/2)+(seqnum*height)-22) : ((height/2)+(seqnum*height)+7);
				var px2start = (locations[i][0] < locations[i][1]) ? (pixel2strand*(locations[i][0]-upstreamstrt)) : (pixel2strand*(locations[i][1]-upstreamstrt));
				var seqboxx = ((0.2*width)+20+settings.padding)+px2start;
				
				var pxofstrand = (locations[i][0] < locations[i][1]) ? (pixel2strand*(locations[i][1]-upstreamstrt)) : (pixel2strand*(locations[i][0]-upstreamstrt));
				var seqboxwidth = pxofstrand - px2start;
				var seqboxheight = 15;
				
				ctx.fillStyle = "#000";
				ctx.fillRect(seqboxx-2, seqboxy-2, seqboxwidth+4, seqboxheight+4);			
				if(typeof window.motifcolors[motif] == 'undefined') {
					do {
						var rndclr = randomcolor();
						var c = rndclr.substring(1);      // strip #
						var rgb = parseInt(c, 16);   // convert rrggbb to decimal
						var r = (rgb >> 16) & 0xff;  // extract red
						var g = (rgb >>  8) & 0xff;  // extract green
						var b = (rgb >>  0) & 0xff;  // extract blue
						var luma = 0.2126 * r + 0.7152 * g + 0.0722 * b; // per ITU-R BT.709
					}
					while(luma < 40 || luma > 200 || rndclr.length < 7);			
					ctx.fillStyle = rndclr;
					window.motifcolors[motif] = rndclr;
				} else {
					ctx.fillStyle = window.motifcolors[motif];
				}				
				ctx.fillRect(seqboxx, seqboxy, seqboxwidth, seqboxheight);
				/* Arrow drawing
				ctx.beginPath();
				var arrowctlpts = [0,3,-10,3,-10,6];
				if (locations[i][2] == "+") {
					ctx.arrow(seqboxx, seqboxy, (seqboxx+seqboxwidth), seqboxy, arrowctlpts);
				} else {
					ctx.arrow((seqboxx+seqboxwidth), seqboxy, seqboxx, seqboxy, arrowctlpts);
				}	*/
				ctx.fill();
				ctx.save();
			}
		}	
		seqnum++;
		ctx.save();
	}
	ctx.save();
	
	// draw legend
	ctx.fillStyle = "#000";
	ctx.fillRect(legendstartx-(settings.padding/2), legendstarty-(settings.padding/2), legendwidth+(settings.padding), legendheight);
	ctx.fillStyle = "#fff";
	ctx.fillRect(legendstartx, legendstarty, legendwidth, legendheight-(settings.padding));
	ctx.fillStyle = "#000";
	ctx.font = settings.legendfont;	
	legendstarty = legendstarty + legdfontsize;
	legendstartx = legendstartx + (settings.padding);
	
	const legendentries = Object.entries(window.motifcolors);
	var legendline = 0;
	for(const [motif, color] of legendentries) {
		if(((legendline+1) % 2) == 0) {
			ctx.fillStyle = color;
			ctx.fillRect(legendstartx+(legendwidth/2), legendstarty+(settings.padding*(legendline-1)+(2*legendline-1))-legdfontsize+2, legdfontsize, legdfontsize);
			ctx.fillStyle = "#000";
			ctx.fillText(motif, legendstartx+legdfontsize+settings.padding+(legendwidth/2), legendstarty+(settings.padding*(legendline-1))+(2*legendline-1));
		} else {
			ctx.fillStyle = color;
			ctx.fillRect(legendstartx, (legendstarty+(settings.padding*legendline)+(2*legendline))-legdfontsize+2, legdfontsize, legdfontsize);
			ctx.fillStyle = "#000";
			ctx.fillText(motif, legendstartx+legdfontsize+settings.padding, legendstarty+(settings.padding*legendline)+(2*legendline));
		}
		legendline++;
	}
}