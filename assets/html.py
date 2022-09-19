# HTML and JS functions

def writeHTML(file, margLikelihoodDict, ipctDict, locusLengthDict):
	with open(file, "w") as outfile:
		outfile.write('''<!DOCTYPE html>
<html lang="en">
	<head>
		<meta charset="utf-8">
		<title>D3js chloroplast gene conflict map</title>
        <script src="https://d3js.org/d3.v4.min.js" charset="utf-8"></script>
        <!--
		<script type="text/javascript" src="genedata.js" charset="utf-8"></script>
		<script type="text/javascript" src="POLPanel.js"></script>
		-->
		<style type="text/css">

			rect:hover {
				fill: orange;
			}

			#tooltip {
				position: absolute;
				width: 200px;
				height: auto;
				padding: 10px;
				background-color: white;
				-webkit-border-radius: 10px;
				-moz-border-radius: 10px;
				border-radius: 10px;
				-webkit-box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
				-moz-box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
				box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
				pointer-events: none;
			}

			#tooltip.hidden {
				display: none;
			}

			#tooltip p {
				margin: 0;
				font-family: sans-serif;
				font-size: 16px;
				line-height: 20px;
			}

		</style>
	</head>
	<body>
		<script type="text/javascript">''')
		outfile.write('''
			var pairs = [''')
		for locuspair in margLikelihoodDict:
			outfile.write('''
			{gene1:"%s", gene2:"%s", symtopo:%s, symtree:%s, apotree1:%s, apotree2:%s},''' % (margLikelihoodDict[locuspair]["gene1"], margLikelihoodDict[locuspair]["gene2"], margLikelihoodDict[locuspair]["concat_brlenUnlinked_ln"], margLikelihoodDict[locuspair]["concat_brlenLinked_ln"], margLikelihoodDict[locuspair]["gene1_ln"], margLikelihoodDict[locuspair]["gene2_ln"]))
		outfile.write('''
			]\n''')
		outfile.write('''
			var genes = [''')
		for gene in ipctDict:
			outfile.write('''
			{name:"%s", chunk:1, ipct:%s, seqlen:%s},''' %(gene, ipctDict[gene], locusLengthDict[gene]))
		outfile.write('''
			]\n''')
		outfile.write('''
			var num_genes = %d''' %(len(ipctDict.keys())))

		outfile.write('''

			function POLPanel(parent, idprefix, t, l, w, h) {
			    this.prefix = idprefix;
			    this.top = t;
			    this.left = l;
			    this.width = w;
			    this.height = h;

			    this.div = parent.append("div")
			        .attr("id", idprefix)
			        .style("position", "absolute")
			        .style("top", t.toString() + "px")
			        .style("left", l.toString() + "px")
			        .style("width", w.toString() + "px")
			        .style("height", h.toString() + "px")
			        .style("vertical-align", "top")
			        .style("text-align", "center");
			    }

			POLPanel.prototype.setup = function() {
			    }

			function POLHelpPanel(parent, idprefix, t, l, w, h) {
			    POLPanel.apply(this, [parent, idprefix, t, l, w, h]);

			    // Immediately hide this.div created in POLPanel constructor
			    this.div.style("display", "none");

			    // Create SVG element
			    this.svg = this.div.append("svg")
			        .attr("id", "helpsvg")
			        .attr("width", this.width)
			        .attr("height", this.height);

			    // define an arrowhead marker
			    // see http://vanseodesign.com/web-design/svg-markers/
			    var defs = this.svg.append("defs");

			    defs.append("marker")
			        .attr("id", "forwardarrow")
			        .attr("markerWidth", 6)
			        .attr("markerHeight", 4)
			        .attr("refX", 0)
			        .attr("refY", 2)
			        .attr("orient", "auto")
			        .attr("markerUnits", "strokeWidth")
			        .append("path")
			        .attr("d", "M0,0 L0,4 L6,2 z")
			        .attr("fill", "white");

			    defs.append("marker")
			        .attr("id", "reversearrow")
			        .attr("markerWidth", 6)
			        .attr("markerHeight", 4)
			        .attr("refX", 6)
			        .attr("refY", 2)
			        .attr("orient", "auto")
			        .attr("markerUnits", "strokeWidth")
			        .append("path")
			        .attr("d", "M6,0 L6,4 L0,2 z")
			        .attr("fill", "white");

			    this.svg.append("rect")
			        .attr("x", 0)
			        .attr("y", 0)
			        .attr("width", this.width)
			        .attr("height", this.height)
			        .attr("fill", d3.color("rgba(0, 0, 0, .6)"))
			        .on("click", function() {
			            d3.select("div#helpbox").style("display", "none");
			        });
			    }

			POLHelpPanel.prototype.infoArrow = function(x0, y0, x1, y1, lw, twoheaded = false) {
			    var arrow = this.svg.append("line")
			        .attr("x1", x0)
			        .attr("y1", y0)
			        .attr("x2", x1)
			        .attr("y2", y1)
			        .attr("stroke", "white")
			        .attr("stroke-width", lw)
			        .attr("marker-end", "url(#forwardarrow)");
			    if (twoheaded)
			        arrow.attr("marker-start", "url(#reversearrow)");
			    }

			POLHelpPanel.prototype.infoText = function(x, y, fontsize, msg) {
			    this.svg.append("text")
			        .attr("x", x)
			        .attr("y", y)
			        .attr("font-family", "Verdana")
			        .attr("font-size", fontsize)
			        .attr("fill", "white")
			        .style("text-anchor", "middle")
			        .text(msg);
			    }

			POLHelpPanel.prototype.infoHTML = function(x, y, fontsize, msg) {
			    this.svg.append("text")
			        .attr("x", x)
			        .attr("y", y)
			        .attr("font-family", "Verdana")
			        .attr("font-size", fontsize)
			        .attr("fill", "white")
			        .style("text-anchor", "middle")
			        .html(msg);
			    }

			function POLCanvasPanel(parent, idprefix, t, l, w, h) {
			    POLPanel.apply(this, [parent, idprefix, t, l, w, h]);

			    // var div = parent.append("div")
			    //     .attr("id", idprefix + "-canvas");

			    // Create SVG element
			    this.svg = this.div.append("svg")
			        .attr("width", this.width)
			        .attr("height", this.height)
			        .on("mousemove", this.mousemoving)
			        .on("mouseout", this.mouseleaving)
			        .on("mousedown", this.dragstarting)
			        .on("mouseup", this.dragstopping);

			    // Draw rect around plot area
			    // this.svg.append("rect")
			    //     .attr("class", "boundingbox")
			    //     .attr("x", 0)
			    //     .attr("y", 0)
			    //     .attr("width", this.width)
			    //     .attr("height", this.height)
			    //     .attr("stroke", "black")
			    //     .attr("fill", "lavender")
			    //     .attr("visibility", (show_bounding_boxes ? "visibility" : "hidden"))
			    //     .style("pointer-events", "none");

			    }

			POLCanvasPanel.prototype.mousemoving = function() {
			    }

			POLCanvasPanel.prototype.mouseleaving = function() {
			    }

			POLCanvasPanel.prototype.dragstarting = function() {
			    }

			POLCanvasPanel.prototype.dragstopping = function() {
			    }

			function POLPlotPanel(parent, idprefix, t, l, w, h, mean, stdev, limits, xlab, ylab) {
			    POLPanel.apply(this, [parent, idprefix, t, l, w, h]);

			    this.dragging_enabled = true;

			    this.xlabel = xlab;
			    this.ylabel = ylab;

			    // Parameters representing mean and standard deviation of the distribution being displayed
			    this.mu = mean;
			    this.sigma = stdev;

			    // Number of segments used to draw curve
			    this.nsegments = 100;

			    // Margins
			    this.left_right_padding = 60;
			    this.top_bottom_padding = 60;

			    // Axis limits
			    this.min_x = limits.xmin;
			    this.max_x = limits.xmax;
			    this.min_y = limits.ymin;
			    this.max_y = limits.ymax;

			    // Create SVG element
			    this.svg = this.div.append("svg")
			        .attr("width", this.width)
			        .attr("height", this.height);

			    // Create scale for X axis (replace with line_scale?)
			    this.xscale = d3.scaleLinear()
			        .domain([this.min_x, this.max_x])
			        .range([this.left_right_padding, this.width - this.left_right_padding]);

			    // Create scale for calculating range bands for x axis
			    this.line_scale = d3.scaleBand()
			        .domain(d3.range(this.nsegments + 1))
			        .range(this.xscale.domain());

			    // Create scale for Y axis
			    this.yscale = d3.scaleLinear()
			        .domain([this.min_y, this.max_y])
			        .range([this.height - this.top_bottom_padding, this.top_bottom_padding]);

			    // Create X axis
			    this.xaxis = d3.axisBottom(this.xscale)
			        .ticks(5)
			        .tickFormat(d3.format(".1f"));

			    // Create Y axis
			    this.yaxis = d3.axisLeft(this.yscale)
			        .ticks(5)
			        .tickFormat(d3.format(".3f"));

			    this.recalcLineData();

			    // Create background rectangle used to capture drag events
			    var bounding_rect = this.svg.append("rect")
			        .attr("x", 0)
			        .attr("y", 0)
			        .attr("width", this.width)
			        .attr("height", this.height)
			        .attr("fill", "white");

			    // Draw hidden rect around plot area that can be made visible for debugging purposes
			    this.svg.append("rect")
			        .attr("class", "boundingbox")
			        .attr("x", 0)
			        .attr("y", 0)
			        .attr("width", this.width)
			        .attr("height", this.height)
			        .attr("stroke", "black")
			        .attr("fill", "lavender")
			        .attr("visibility", "hidden")
			        .style("pointer-events", "none");

			    // Add X axis to svg
			    this.svg.append("g")
			        .attr("id", this.prefix + "-xaxis")
			        .attr("class", "axis")
			        .attr("transform", "translate(0," + (this.height - this.top_bottom_padding) + ")")
			        .call(this.xaxis);

			    // Add Y axis to svg
			    this.svg.append("g")
			        .attr("id", this.prefix + "-yaxis")
			        .attr("class", "axis")
			        .attr("transform", "translate(" + this.left_right_padding + ",0)")
			        .call(this.yaxis);

			    // Style the axes
			    this.svg.selectAll('.axis line, .axis path')
			        .style({'stroke': 'black', 'fill': 'none', 'stroke-width': '2px', 'shape-rendering': 'crispEdges'});
			    this.svg.selectAll('.axis text')
			        .style({'font-family': 'Verdana', 'font-size': '28px'});
			    this.svg.selectAll('#xaxis.axis text')
			        .attr("transform", "translate(0, 10)");

			    // Create x axis label
			    if (this.xlabel !== null) {
			        this.svg.append("text")
			            .attr("id", this.prefix + "-xlab")
			            .attr("x", this.width/2)
			            .attr("y", this.height - 20)
			            .style("text-anchor", "middle")
			            .style("font-family", "Verdana")
			            .style("font-size", "8pt")
			            .text(this.xlabel);
			        }

			    // Create y axis label
			    if (this.ylabel !== null) {
			        this.svg.append("text")
			            .attr("id", this.prefix + "-ylab")
			            .attr("x", 10)
			            .attr("y", this.height/2)
			            .attr("transform", "rotate(-90 10," + this.height/2 + ")")
			            .style("text-anchor", "middle")
			            .style("font-family", "Verdana")
			            .style("font-size", "8pt")
			            .text(this.ylabel);
			        }

			    var self = this;

			    // Create line generator function
			    this.lineFunc = d3.line()
			        .x(function(d) {return self.xscale(d.x);})
			        .y(function(d) {return self.yscale(d.y);});

			    // Create line representing probability density
			    this.density = this.svg.append("path")
			        .attr("d", this.lineFunc(this.linedata))
			        .attr("fill", "none")
			        .attr("stroke", "blue")
			        .attr("stroke-width", 2)
			        .style("pointer-events", "none");   // don't want line intercepting drag events

			    // Create scales for choosing scaling factors for sigma based on drag extent
			    // The sigma_pos_scale is used if user drags downward
			    // The sigma_pos_scale is used if user drags upward
			    var sigma_pos_scale = d3.scaleLinear()
			        .domain([0,this.height])
			        .range([0.0,limits.dymouse[1]]);
			    var sigma_neg_scale = d3.scaleLinear()
			        .domain([0,-this.height])
			        .range([0.0,limits.dymouse[0]]);

			    // Create drag behavior
			    y_at_drag_start = null;
			    sigma_at_drag_start = null;
			    drag = d3.drag()
			        .on("start", function(d) {
			            if (self.dragging_enabled) {
			                y_at_drag_start = d3.event.y;
			                sigma_at_drag_start = self.sigma;
			                d3.event.sourceEvent.stopPropagation();
			                d3.select(this).classed("dragging", true);
			                }
			        })
			        .on("drag", function(d) {
			            if (self.dragging_enabled) {
			                // Move mu by the amount corresponding to the x-component of the drag
			                self.mu = self.xscale.invert(self.xscale(self.mu) + d3.event.dx);
			                if (self.mu < 0.0001)
			                    self.mu = 0.0001;

			                // Adjust sigma based on extent of drag in vertical direction
			                var dy = d3.event.y - y_at_drag_start;
			                if (dy < 0) {
			                    self.sigma = sigma_at_drag_start*Math.exp(sigma_neg_scale(dy));
			                } else {
			                    self.sigma = sigma_at_drag_start*Math.exp(sigma_pos_scale(dy));
			                }

			                self.refreshPlot();
			                self.duringdrag();
			                }
			        })
			        .on("end", function(d) {
			            if (self.dragging_enabled) {
			                d3.select(this).classed("dragging", false);

			                // Recalculate linedata
			                self.recalcLineData();
			                self.density.transition()
			                    .duration(500)
			                    .attr("d", self.lineFunc(self.linedata));

			                self.afterdrag();
			                }
			        });

			    bounding_rect.call(drag);

			    }

			POLPlotPanel.prototype.enableDragging = function() {
			    this.dragging_enabled = true;
			    console.log("dragging enabled");
			    }

			POLPlotPanel.prototype.disableDragging = function() {
			    this.dragging_enabled = false;
			    console.log("dragging disabled");
			    }

			POLPlotPanel.prototype.duringdrag = function() {
			    }

			POLPlotPanel.prototype.afterdrag = function() {
			    }

			POLPlotPanel.prototype.refreshPlot = function() {
			    // Recalculate linedata
			    this.recalcLineData();

			    // Recalculate path
			    this.density.attr("d", this.lineFunc(this.linedata));
			    }

			POLPlotPanel.prototype.calcLogDensity = function(x) {
			    // Exponential with mean mu
			    var lambda = 1.0/this.param2;
			    var logy = Math.log(lambda) - lambda*x;
			    return logy;
			    }

			POLPlotPanel.prototype.recalcLineData = function() {
			    this.recalcParams();    // calculate param1 and param2 from mu and sigma
			    this.linedata = [];
			    for (var i = 1; i < this.nsegments + 1; i++) {
			        // note that we skip i=0 to avoid calculating density at x=0.0 (which may be infinity)
			        var x = this.line_scale(i);
			        var logy = this.calcLogDensity(x);
			        var y = Math.exp(logy);
			        this.linedata.push({'x':x, 'y':y});
			        //console.log(this.prefix + ": i = " + i + ", x = " + x + ", y = " + y + ", logy = " + logy);
			        }
			    }

			POLPlotPanel.prototype.recalcParams = function() {
			    this.param1 = 1.0;
			    this.param2 = this.mu;
			    //console.log("recalcParams: mu = " + this.mu + ", sigma = " + this.sigma + ", shape = " + this.param1 + ", scale = " + this.param2);
			    }

			var addStringDropdown = function(panel, label, choices, default_choice, onfunc) {
			    var control_div = panel.append("div").append("div")
			        .attr("class", "control");
			    var select = control_div.append("select")
			        .on("change", onfunc);
			    select.selectAll("option")
			        .data(choices)
			        .enter()
			        .append("option")
			        .property("selected", function(d,i) {return (i == default_choice ? true : false);})
			        .text(function(d) {return d;});
			    control_div.append("label")
			        .style("font-family", "Verdana")
			        .style("font-size", "10pt")
			        .html("&nbsp;" + label);
			    //console.log("selectedIndex for \"" + label + "\" dropdown = " + select.property("selectedIndex"));
			    }

			var addIntDropdown = function(panel, label, choices, default_choice, onfunc) {
			    var control_div = panel.append("div").append("div")
			        .attr("class", "control");
			    var select = control_div.append("select")
			        .on("change", onfunc);
			    select.selectAll("option")
			        .data(choices)
			        .enter()
			        .append("option")
			        .property("selected", function(d,i) {return (i == default_choice ? true : false);})
			        .text(function(d) {return d.toFixed(0);});
			    control_div.append("label")
			        .style("font-family", "Verdana")
			        .style("font-size", "10pt")
			        .html("&nbsp;" + label);
			    //console.log("selectedIndex for \"" + label + "\" dropdown = " + select.property("selectedIndex"));
			    }

			var addCheckbox = function(panel, label, onfunc) {
			    var control_div = panel.append("div").append("div")
			        .attr("class", "control");
			    control_div.append("input")
			        .attr("type", "checkbox")
			        .on("change", onfunc);
			    control_div.append("label")
			        .append("label")
			        .html("&nbsp;" + label);
			    }

			var addRadioButtons = function(panel, name, values, default_value, label, onfunc) {
			    var control_div = panel.append("div").append("div")
			        .attr("class", "control");
			    control_div.append("label")
			        .append("label")
			        .html(label + "&nbsp;");
			    for (v in values) {
			        control_div.append("input")
			            .attr("type", "radio")
			            .attr("name", name)
			            .attr("value", values[v])
			            .property("checked", (values[v] == default_value ? true : false))
			            .on("change", onfunc);
			        control_div.append("label")
			            .html(values[v]);
			        }
			    }

			var addButton = function(panel, label, onfunc) {
			    var control_div = panel.append("div")
			        .attr("class", "control");
			    control_div.append("input")
			        .attr("type", "button")
			        .attr("value", label)
			        .on("click", onfunc);
			    }


''')

		outfile.write('''
		    var color_bar_height = 30;

            function inherit(p) {
                if (p == null)
                    throw TypeError();
                if (Object.create)
                    return Object.create(p);
                var t = typeof p;
                if (t !== "object" && t !== "function")
                    throw TypeError();
                function f() {};
                f.prototype = p;
                return new f();
                }

            // modified from https://krazydad.com/tutorials/makecolors.php
            function makeColorGradient(frequency1, frequency2, frequency3, phase1, phase2, phase3, center, width, len) {
                if (center == undefined)
                    center = 128;
                if (width == undefined)
                    width = 127;
                if (len == undefined)
                    len = 50;

                var colors = [];
                for (var i = 0; i < len; ++i) {
                    var r = Math.sin(frequency1*i + phase1) * width + center;
                    var g = Math.sin(frequency2*i + phase2) * width + center;
                    var b = Math.sin(frequency3*i + phase3) * width + center;
                    var color = d3.rgb(r, g, b);
                    colors.push(color);
                    }
                return colors;
                }

            var ncolors = 10;
            var colordata = makeColorGradient(-Math.PI/ncolors, 0, Math.PI/ncolors, 3.*Math.PI/2, 3.*Math.PI/2, Math.PI/2, 128, 127, ncolors);
            var color_scale = d3.scaleOrdinal(colordata).domain(d3.range(ncolors));
            //var color_scale = d3.scaleOrdinal(d3.schemeCategory20b);

            function GeneCirclePanel(parent, prefix, top, left, width, height) {
                POLCanvasPanel.apply(this, arguments);
                }
            GeneCirclePanel.prototype = inherit(POLCanvasPanel.prototype);
            GeneCirclePanel.constructor = GeneCirclePanel;

            GeneCirclePanel.prototype.setup = function() {
                this.plot_width = this.width;
                this.plot_height = this.height - color_bar_height;

                this.outer_radius = this.plot_width/2 - 50;
                this.inner_radius = this.outer_radius - 50;
                this.label_radius = this.outer_radius + 10;
                this.tension = 1.0;
                this.drag_start = null;

                var minimum_value = 0;
                var maximum_value = 100;

                this.color_genes_by_info = true;

                var formatAsPercentage = d3.format(".1");

                this.createRadialLine = function(tension) {
                    return d3.radialLine()
                        .angle(function(d) {return d.angle;})
                        .radius(function(d) {return d.radius;})
                        .curve(d3.curveBundle.beta(tension))
                    }

                rl = this.createRadialLine(1.0);

                var labelXY = function(angleRadians, radius) {
                    var dx = Math.cos(angleRadians)*radius;
                    var dy = Math.sin(angleRadians)*radius;
                    return [dx,dy];
                    }

                var separate_better_than_common_topo = function(d, tol) {
                    return (d.apotree1 + d.apotree2) - d.symtopo > tol;
                    }

                var common_topo_better_than_separate = function(d, tol) {
                    return d.symtopo - (d.apotree1 + d.apotree2) > tol;
                    }

                var separate_better_than_common_tree = function(d, tol) {
                    return (d.apotree1 + d.apotree2) - d.symtree > tol;
                    }

                var common_tree_better_than_separate = function(d, tol) {
                    return d.symtree - (d.apotree1 + d.apotree2) > tol;
                    }

                var common_tree_better_than_common_topo = function(d, tol) {
                    return d.symtree - d.symtopo > tol;
                    }

                var common_topo_better_than_common_tree = function(d, tol) {
                    return d.symtree - d.symtopo > tol;
                    }

                var show_connection = function(d) {
                    var tol = 1.0;

                    return separate_better_than_common_topo(d, tol);
                    //return common_topo_better_than_separate(d, tol);

                    //return separate_better_than_common_tree(d, tol);
                    //return common_tree_better_than_separate(d, tol);

                    //return common_tree_better_than_common_topo(d, tol);
                    //return common_topo_better_than_common_tree(d, tol);
                    }

                // eps is the space (subtending angle) between chunks
                var eps = .02;

                // calculate total length of all genes
                var nsites_total = d3.sum(genes, function(d) {return d.seqlen;});

                // count number of chunks (genes in same chunk have never been separated)
                var nchunks = d3.max(genes, function(d) {return d.chunk;});

                // total_angle is fraction of circle devoted to chunks, excluding gaps between chunks
                var total_angle = 2.147*Math.PI - eps*nchunks;

                // store information for arcs representing genes in arcinfo

                var genes_in_order = [];
								for (var k = 0; k < num_genes; k++) {
									genes_in_order.push(genes[k].name)
								}
								// new Array(function(d) {return d.name;});
                var arcinfo = [];
                var cum = 0;
                for (var g = 0; g < genes_in_order.length; g++) {
                    var i = genes.findIndex(function(d) {return d.name == genes_in_order[g];});
                    var ipct = parseInt(genes[i].ipct);
                    var seqlen = parseInt(genes[i].seqlen);
                    var startfrac     = 1.0*cum/nsites_total;
                    var endfrac       = 1.0*(cum + seqlen)/nsites_total;
                    var start_angle   = total_angle*startfrac;
                    var end_angle     = total_angle*endfrac;
                    var arcargs = {
                        gene:        genes[i].name,
                        seqlen:      seqlen,
                        ipct:        ipct,
                        chunk:       parseInt(genes[i].chunk),
                        color_index: 0,
                        innerRadius: this.inner_radius,
                        outerRadius: this.outer_radius,
                        startAngle:  start_angle,
                        endAngle:    end_angle,
                        leftmost:    false,
                        rightmost:   false,
                        linkages:    []             // filled in by loop over pairs below
                        };
                    arcinfo.push(arcargs);
                    cum += seqlen;
                    }

                // fix up chunk boundaries to create small gap between chunks
                var gene_weight = {};   // records number of other genes to which key gene is connected
                var gene_angle = [];
                var cum_eps = 0.0;
                for (var i in arcinfo) {
                    var j = (i > 0 ? i : arcinfo.length) - 1;
                    if (arcinfo[i].chunk != arcinfo[j].chunk) {
                        cum_eps += eps;
                        arcinfo[i].leftmost = true;
                        arcinfo[j].rightmost = true;
                        }
                    arcinfo[i].startAngle += cum_eps;
                    arcinfo[i].endAngle   += cum_eps;

                    // calculate angle and position of gene label
                    var average_angle = (arcinfo[i].startAngle + arcinfo[i].endAngle)/2;
                    gene_angle.push({gene:arcinfo[i].gene, angle:average_angle, radius:this.inner_radius});
                    arcinfo[i].labelpos = labelXY(average_angle - Math.PI/2.0, this.label_radius);
                    arcinfo[i].labelangle = 180.0*(average_angle - Math.PI/2.0)/Math.PI;
                    arcinfo[i].labelflip = (arcinfo[i].labelangle > 90 && arcinfo[i].labelangle < 270);

                    // initialize gene weight
                    gene_weight[arcinfo[i].gene] = 0;
                    }

                // store information about pairs in pairinfo
                var pairinfo = [];
                for (i = 0; i < pairs.length; i++) {
                    var index_gene1 = gene_angle.findIndex(function(d) {return d.gene == pairs[i].gene1;});
                    var index_gene2 = gene_angle.findIndex(function(d) {return d.gene == pairs[i].gene2;});
                    var pinfo = {
                        gene1: pairs[i].gene1,
                        gene2: pairs[i].gene2,
                        index1: index_gene1,
                        index2: index_gene2,
                        chunk1: genes[genes.findIndex(function(d) {return d.name == pairs[i].gene1;})].chunk,
                        chunk2: genes[genes.findIndex(function(d) {return d.name == pairs[i].gene2;})].chunk,
                        angle1: gene_angle[index_gene1].angle,
                        angle2: gene_angle[index_gene2].angle,
                        radius1: gene_angle[index_gene1].radius,
                        radius2: gene_angle[index_gene2].radius,
                        symtopo: parseFloat(pairs[i].symtopo),
                        symtree: parseFloat(pairs[i].symtree),
                        apotree1: parseFloat(pairs[i].apotree1),
                        apotree2: parseFloat(pairs[i].apotree2)
                        };
                    pairinfo.push(pinfo);

                    if (show_connection(pinfo)) {
                        // increment weights for each gene in pair
                        gene_weight[pairinfo[i].gene1]++;
                        gene_weight[pairinfo[i].gene2]++;

                        // add gene linkages to relevant arcinfo elements
                        index_gene1 = arcinfo.findIndex(function(d) {return d.gene == pairs[i].gene1;});
                        index_gene2 = arcinfo.findIndex(function(d) {return d.gene == pairs[i].gene2;});
                        arcinfo[index_gene1].linkages.push(pairs[i].gene2);
                        arcinfo[index_gene2].linkages.push(pairs[i].gene1);
                        }
                    }

                // find minimum and maximum gene weights (gene weight is number of connections a gene has to other genes)
                var min_gene_weight = d3.min(d3.values(gene_weight));
                var max_gene_weight = d3.max(d3.values(gene_weight));
                var gene_weight_span = max_gene_weight - min_gene_weight;
                console.log("min_gene_weight = " + min_gene_weight);
                console.log("max_gene_weight = " + max_gene_weight);

                // find minimum and maximum information content
                var min_gene_info = d3.min(arcinfo, function(d) {return d.ipct;});
                var max_gene_info = d3.max(arcinfo, function(d) {return d.ipct;});
                var info_span = max_gene_info - min_gene_info;
                console.log("min_gene_info = " + min_gene_info);
                console.log("max_gene_info = " + max_gene_info);

                // fix up arcinfo to specify color_scale indices
                for (var i in arcinfo) {
                    if (this.color_genes_by_info) {
                        var I = arcinfo[i].ipct;
                        var c = parseInt(1.0*ncolors*(I - min_gene_info)/info_span);
                        }
                    else {
                        var w = gene_weight[arcinfo[i].gene]
                        var c = parseInt(1.0*ncolors*(w - min_gene_weight)/gene_weight_span);
                    }
                    arcinfo[i].color_index = (c < ncolors ? c : ncolors - 1);
                    //console.log("gene = " + arcinfo[i].gene + ", w = " + w + ", color_index = " + arcinfo[i].color_index);
                    }

                var self = this;

                this.svg.append("rect")
                    .attr("id", "plotarea")
                    .attr("x", 0)
                    .attr("y", 0)
                    .attr("width", this.plot_width)
                    .attr("height", this.plot_height)
                    .style("fill", "lavender")
                    .attr("visibility", "visible");

                this.svg.append("circle")
                    .attr("cx", this.plot_width/2)
                    .attr("cy", this.plot_height/2)
                    .attr("r", 2)
                    .style("fill", "black");

                this.svg.selectAll("rect.colorspot")
                    .data(colordata)
                    .enter()
                    .append("rect")
                    .attr("class", "colorspot")
                    .attr("x", function(d,i) {return i*self.plot_width/ncolors;})
                    .attr("y", this.plot_height)
                    .attr("width", function(d,i) {return self.plot_width/ncolors;})
                    .attr("height", color_bar_height)
                    .attr("visibility", "visible")
                    .style("fill", function(d) {return d;})
                    .style("stroke", "white");

                this.svg.selectAll("text.colorspot")
                    .data(colordata)
                    .enter()
                    .append("text")
                    .attr("class", "colorspot")
                    .attr("x", function(d,i) {return i*self.plot_width/ncolors + 0.5*self.plot_width/ncolors;})
                    .attr("y", this.plot_height + 0.7*color_bar_height)
                    .attr("width", function(d,i) {return self.plot_width/ncolors;})
                    .attr("height", color_bar_height)
                    .attr("visibility", "visible")
                    .style("font-family", "Verdana")
                    .style("font-size", "10pt")
                    .style("fill", "white")
                    .style("text-anchor", "middle")
                    .text(function(d,i) {
                        if (self.color_genes_by_info) {
                            var lower = parseInt(1.0*i*info_span/ncolors + min_gene_info);
                            var upper = parseInt(1.0*(i+1)*info_span/ncolors + min_gene_info);
                            return lower + "-" + upper;
                            }
                        else {
                            var lower = parseInt(1.0*i*gene_weight_span/ncolors + min_gene_weight);
                            var upper = parseInt(1.0*(i+1)*gene_weight_span/ncolors + min_gene_weight);
                            return lower + "-" + upper;
                            }
                        });

                // create arcs representing genes
                var gene_arcs = this.svg.append("g")
                    .attr("id", "gene_arcs");
                gene_arcs.selectAll("path")
                    .data(arcinfo)
                    .enter()
                    .append("path")
                    .attr("class", function(d) {return "gene-" + d.gene;})
                    .attr("d", function(d) {
                        return d3.arc()(d);
                        })
                    //.attr("fill", function(d,i) {return color_scale(i)})
                    .attr("fill", function(d) {return color_scale(d.color_index)})
                    .attr("transform", "translate(" + (this.plot_width/2) + "," + (this.plot_height/2) + ")")
                    .on("mouseover", function(d) {
                        //Update the tooltip position and value
                        d3.select("div#tooltip")
                            .style("left", parseInt(d.labelpos[0] + self.plot_width/2) + "px")
                            .style("top",  parseInt(d.labelpos[1] + self.plot_height/2) + "px");
                        d3.select("span#genename")
                            .html(d.gene);
                        d3.select("span#geneinfo")
                            .html("Ipct = " + d.ipct + "<br/>seqlen = " + d.seqlen);

                        //Show the tooltip
                        d3.select("div#tooltip").classed("hidden", false);

                        d3.select(this)
                            .attr("fill", "orange");
                            /*.attr("d", d3.arc()({
                                innerRadius: d.innerRadius - 5,
                                outerRadius: d.outerRadius + 5,
                                startAngle:  d.startAngle,
                                endAngle:    d.endAngle}));*/
                        gene_links.selectAll("path")
                            .attr("visibility", "hidden");
                        gene_links.selectAll("path." + d.gene)
                            .filter(function(d) {return d.chunk1 == d.chunk2;})
                            .attr("visibility", "visible")
                            .style("stroke-width", "3px")
                            .style("stroke", "red");
                        gene_links.selectAll("path." + d.gene)
                            .filter(function(d) {return d.chunk1 != d.chunk2;})
                            .attr("visibility", "visible")
                            .style("stroke-width", "3px")
                            .style("stroke", "blue");
                        })
                    .on("mouseout", function(d,i) {
                        //Hide the tooltip
                        d3.select("div#tooltip").classed("hidden", true);

                        d3.select(this)
                            //.attr("fill", color_scale(i))
                            .attr("fill", function(d) {return color_scale(d.color_index)})
                            .attr("d", d3.arc()(d));
                        gene_links.selectAll("path")
                            .attr("visibility", "visible")
                            .style("stroke-width", "1px")
                            .style("stroke", "gray");
                        });

                // create lines separating genes
                gene_arcs.selectAll("line.leftsep")
                    .data(arcinfo)
                    .enter()
                    .append("line")
                    .attr("class", "leftsep")
                    .attr("x1", function(d) {return d.innerRadius*Math.cos(d.startAngle - Math.PI/2);})
                    .attr("y1", function(d) {return d.innerRadius*Math.sin(d.startAngle - Math.PI/2);})
                    .attr("x2", function(d) {return d.outerRadius*Math.cos(d.startAngle - Math.PI/2);})
                    .attr("y2", function(d) {return d.outerRadius*Math.sin(d.startAngle - Math.PI/2);})
                    .attr("visibility", "visible")
                    .style("stroke", "black")
                    .style("stroke-width", "1px")
                    .attr("transform", "translate(" + (this.plot_width/2) + "," + (this.plot_height/2) + ")")

                gene_arcs.selectAll("line.rightsep")
                    .data(arcinfo)
                    .enter()
                    .append("line")
                    .attr("filter", function(d) {return d.rightmost;})
                    .attr("class", "rightsep")
                    .attr("x1", function(d) {return d.innerRadius*Math.cos(d.endAngle - Math.PI/2);})
                    .attr("y1", function(d) {return d.innerRadius*Math.sin(d.endAngle - Math.PI/2);})
                    .attr("x2", function(d) {return d.outerRadius*Math.cos(d.endAngle - Math.PI/2);})
                    .attr("y2", function(d) {return d.outerRadius*Math.sin(d.endAngle - Math.PI/2);})
                    .attr("visibility", "visible")
                    .style("stroke", "black")
                    .style("stroke-width", "1px")
                    .attr("transform", "translate(" + (this.plot_width/2) + "," + (this.plot_height/2) + ")");

                // create labels for genes
                var gene_labels = this.svg.append("g")
                    .attr("id", "gene_labels");
                gene_labels.selectAll("text")
                    .data(arcinfo)
                    .enter()
                    .append("text")
                    .attr("x", function(d) {return d.labelpos[0];})
                    .attr("y", function(d) {return d.labelpos[1];})
                    .attr("transform", function(d) {
                        var rotate_angle = d.labelflip ? (d.labelangle-180) : d.labelangle;
                        var rotate_about_x = d.labelpos[0] + self.plot_width/2;
                        var rotate_about_y = d.labelpos[1] + self.plot_height/2;
                        return "rotate(" + rotate_angle + "," + rotate_about_x + "," + rotate_about_y + ") translate(" + (self.plot_width/2) + "," + (self.plot_height/2) + ")";
                        })
                    .attr("fill", "black")
                    .style("font-family", "verdana")
                    .style("font-size", "10pt")
                    .attr("text-anchor", function(d) {return (d.labelflip ? "end" : "start");})
                    .text(function(d) {return d.gene;});

                // create curved path linking genes
                var gene_links = this.svg.append("g")
                    .attr("id", "gene_links");
                gene_links.selectAll("path")
                    .data(pairinfo)
                    .enter()
                    .filter(function(d) {return show_connection(d);})
                    .append("path")
                    .attr("class", function(d) {return "lnks " + d.gene1 + " " + d.gene2;})
                    .attr("d", function(d) {
                        var line_data = [{angle:d.angle1,radius:d.radius1},{angle:0, radius:0},{angle:d.angle2,radius:d.radius2}];
                        return rl(line_data);
                        })
                    .style("fill", "none")
                    .style("stroke", "gray")
                    .style("stroke-width", "1px")
                    .attr("transform", "translate(" + (this.plot_width/2) + "," + (this.plot_height/2) + ")");
                }

            GeneCirclePanel.prototype.setDragStartPoint = function(mousepos) {
                this.drag_start = mousepos;
                }

            GeneCirclePanel.prototype.clearDragStartPoint = function() {
                this.drag_start = null;
                }

            GeneCirclePanel.prototype.dragstarting = function() {
                gene_circle.setDragStartPoint(d3.mouse(this));
                }

            GeneCirclePanel.prototype.dragstopping = function() {
                gene_circle.clearDragStartPoint();
                }

            GeneCirclePanel.prototype.mousemoving = function() {
                if (gene_circle.drag_start) {
                    var y0 = gene_circle.drag_start[1];
                    var coords = d3.mouse(this);
                    var y = coords[1];
                    var frac = (y0 - y)/gene_circle.plot_height;
                    var prev_tension = gene_circle.tension;
                    gene_circle.tension += frac*(y > y0 ? gene_circle.tension : 1.0 - gene_circle.tension);

                    var rl = gene_circle.createRadialLine(gene_circle.tension);
                    d3.selectAll("path.lnks")
                        .attr("d", function(d) {
                            var line_data = [{angle:d.angle1,radius:d.radius1},{angle:0, radius:0},{angle:d.angle2,radius:d.radius2}];
                            return rl(line_data);
                            });

                    this.drag_start = coords;
                    }
                }

            GeneCirclePanel.prototype.mouseleaving = function() {
                //console.log("mouse leaving");
                }

            var container_div = d3.select("body");
            var gene_circle = new GeneCirclePanel(container_div, "genecircle", 0, 0, 800, 800 + color_bar_height);
            gene_circle.setup();

		</script>
		<div id="tooltip" class="hidden">
			<p><strong><span id="genename">xxx</span></strong></p>
			<p><span id="geneinfo">xxx</span></p>
		</div>
	</body>
</html>


''')
