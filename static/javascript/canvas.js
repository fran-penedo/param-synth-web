// from http://stackoverflow.com/questions/610406/javascript-equivalent-to-printf-string-format
if (!String.prototype.format) {
  String.prototype.format = function() {
    var args = arguments;
    return this.replace(/{(\d+)}/g, function(match, number) { 
      return typeof args[number] != 'undefined'
        ? args[number]
        : match
      ;
    });
  };
}

function Graph(config) {
    // user defined properties
    var qcanvas = $("#" + config.canvasId);
    this.canvas = qcanvas.get(0);
    this.minX = config.minX;
    this.minY = config.minY;
    this.maxX = config.maxX;
    this.maxY = config.maxY;
    this.unitsPerTick = config.unitsPerTick;
    qcanvas.click(this, handleClick);
    qcanvas.dblclick(this, handleDblClick);
    this.poly = [];
    this.poly_closed = false;
    this.cur_constr = [];
    this.constrs = [];
    if (config.main) {
        this.par_canvas = new Graph({
            canvasId: 'par_canvas',
            minX: -10,
            minY: -10,
            maxX: 10,
            maxY: 10,
            unitsPerTick: 1
        });
        this.par_canvas.handler = Handlers["poly"];
    }

    // constants
    this.axisColor = '#aaa';
    this.font = '8pt Calibri';
    this.tickSize = 10;

    // relationships
    this.context = this.canvas.getContext('2d');
    this.rangeX = this.maxX - this.minX;
    this.rangeY = this.maxY - this.minY;
    this.unitX = this.canvas.width / this.rangeX;
    this.unitY = this.canvas.height / this.rangeY;
    this.centerY = Math.round(Math.abs(this.minY / this.rangeY) * this.canvas.height);
    this.centerX = Math.round(Math.abs(this.minX / this.rangeX) * this.canvas.width);
    this.iteration = (this.maxX - this.minX) / 1000;
    this.scaleX = this.canvas.width / this.rangeX;
    this.scaleY = this.canvas.height / this.rangeY;


    // draw x and y axis
    this.erase();
}

Graph.prototype.drawXAxis = function() {
    var context = this.context;
    context.save();
    context.beginPath();
    context.moveTo(0, this.centerY);
    context.lineTo(this.canvas.width, this.centerY);
    context.strokeStyle = this.axisColor;
    context.lineWidth = 2;
    context.stroke();

    // draw tick marks
    var xPosIncrement = this.unitsPerTick * this.unitX;
    var xPos, unit;
    context.font = this.font;
    context.textAlign = 'center';
    context.textBaseline = 'top';

    // draw left tick marks
    xPos = this.centerX - xPosIncrement;
    unit = -1 * this.unitsPerTick;
    while(xPos > 0) {
        context.moveTo(xPos, this.centerY - this.tickSize / 2);
        context.lineTo(xPos, this.centerY + this.tickSize / 2);
        context.stroke();
        //context.fillText(unit, xPos, this.centerY + this.tickSize / 2 + 3);
        unit -= this.unitsPerTick;
        xPos = Math.round(xPos - xPosIncrement);
    }

    // draw right tick marks
    xPos = this.centerX + xPosIncrement;
    unit = this.unitsPerTick;
    while(xPos < this.canvas.width) {
        context.moveTo(xPos, this.centerY - this.tickSize / 2);
        context.lineTo(xPos, this.centerY + this.tickSize / 2);
        context.stroke();
        if (unit == 1) {
            context.fillText(unit, xPos, this.centerY + this.tickSize / 2 + 3);
        }
        unit += this.unitsPerTick;
        xPos = Math.round(xPos + xPosIncrement);
    }
    context.restore();
};

Graph.prototype.drawYAxis = function() {
    var context = this.context;
    context.save();
    context.beginPath();
    context.moveTo(this.centerX, 0);
    context.lineTo(this.centerX, this.canvas.height);
    context.strokeStyle = this.axisColor;
    context.lineWidth = 2;
    context.stroke();

    // draw tick marks
    var yPosIncrement = this.unitsPerTick * this.unitY;
    var yPos, unit;
    context.font = this.font;
    context.textAlign = 'right';
    context.textBaseline = 'middle';

    // draw top tick marks
    yPos = this.centerY - yPosIncrement;
    unit = this.unitsPerTick;
    while(yPos > 0) {
        context.moveTo(this.centerX - this.tickSize / 2, yPos);
        context.lineTo(this.centerX + this.tickSize / 2, yPos);
        context.stroke();
        if (unit == 1) {
            context.fillText(unit, this.centerX - this.tickSize / 2 - 3, yPos);
        }
        unit += this.unitsPerTick;
        yPos = Math.round(yPos - yPosIncrement);
    }

    // draw bottom tick marks
    yPos = this.centerY + yPosIncrement;
    unit = -1 * this.unitsPerTick;
    while(yPos < this.canvas.height) {
        context.moveTo(this.centerX - this.tickSize / 2, yPos);
        context.lineTo(this.centerX + this.tickSize / 2, yPos);
        context.stroke();
        //context.fillText(unit, this.centerX - this.tickSize / 2 - 3, yPos);
        unit -= this.unitsPerTick;
        yPos = Math.round(yPos + yPosIncrement);
    }
    context.restore();
};

Graph.prototype.drawEquation = function(equation, color, thickness) {
    var context = this.context;
    context.save();
    context.save();
    this.transformContext();

    context.beginPath();
    context.moveTo(this.minX, equation(this.minX));

    for(var x = this.minX + this.iteration; x <= this.maxX; x += this.iteration) {
        context.lineTo(x, equation(x));
    }

    context.restore();
    context.lineJoin = 'round';
    context.lineWidth = thickness;
    context.strokeStyle = color;
    context.stroke();
    context.restore();
};

Graph.prototype.transformContext = function() {
    var context = this.context;

    // move context to center of canvas
    this.context.translate(this.centerX, this.centerY);

    /*
        * stretch grid to fit the canvas window, and
        * invert the y scale so that that increments
        * as you move upwards
        */
    context.scale(this.scaleX, -this.scaleY);
};

Graph.prototype.toScale = function(x, y) {
    return [(x - this.centerX) / this.scaleX, - (y - this.centerY) / this.scaleY]
};

Graph.prototype.fromScale = function(x, y) {
    return [(x  * this.scaleX + this.centerX), (- y * this.scaleY + this.centerY)]
};


Graph.prototype.drawPoly = function(poly, fill) {
    fill = typeof fill !== 'undefined' ? fill : false;
    var p = this.poly;
    this.poly = poly;
    this.poly_closed = true;
    this.renderPoly(fill);
    this.poly = p;
};

Graph.prototype.erase = function() {
    this.context.clearRect(0, 0, canvas.width, canvas.height);
    this.drawXAxis();
    this.drawYAxis();
};

Graph.prototype.renderPoly = function(fill) {
    fill = typeof fill !== 'undefined' ? fill : false;
    var context = this.context;
    context.save();
    this.transformContext();
    context.lineWidth = 0.1;

    if (this.poly != undefined && this.poly.length > 1) {
        context.beginPath();
        context.moveTo(this.poly[0][0], this.poly[0][1]);
        for (var i = 1; i < this.poly.length; i++) {
            context.lineTo(this.poly[i][0], this.poly[i][1]);
        }
        if (this.poly_closed) {
            context.closePath();
        }
        if (fill) {
            this.context.fillStyle = 'rgba(0, 255, 0, 0.2)'
            context.fill();
        } else {
            context.stroke();
        }
    }
    context.restore();
};

Graph.prototype.renderCurConstr = function() {
    var context = this.context;
    context.save();
    this.transformContext();
    context.lineWidth = 0.1;

    if (this.cur_constr.length > 1) {
        var a = this.cur_constr[0],
            b = this.cur_constr[1],
            slope = (b[1] - a[1]) / (b[0] - a[0]),
            yinterl = a[1] + slope * (this.minX - a[0]),
            yinterr = a[1] + slope * (this.maxX - a[0]),
            xinterd = (this.minY - a[1]) / slope + a[0],
            xinteru = (this.maxY - a[1]) / slope + a[0];

        context.beginPath();

        if (yinterl < this.minY) {
            context.moveTo(xinterd, this.minY);
        } else if (yinterl > this.maxY){
            context.moveTo(xinteru, this.maxY);
        } else {
            context.moveTo(this.minX, yinterl);
        }

        if (yinterr < this.minY) {
            context.lineTo(xinterd, this.minY);
        } else if (yinterr > this.maxY){
            context.lineTo(xinteru, this.maxY);
        } else {
            context.lineTo(this.maxX, yinterr);
        }
        context.stroke();
    }

    context.restore();
};

Graph.prototype.renderConstrs = function() {
    for (var i = 0; i < this.constrs.length; i++) {
        this.cur_constr = this.constrs[i];
        this.renderCurConstr();
    }
};

Graph.prototype.drawPartitionNames = function() {
    var context = this.context;
    context.save();

    for (var i = 0; i < this.partition.length; i++) {
        var set = this.partition[i];
        var center = set["centroid"];
        var tcenter = this.fromScale(center[0], center[1]);
        context.fillText(set["name"], tcenter[0], tcenter[1]);
    }

    context.restore();
};

Graph.prototype.click = function(x, y) {
    this.handler.click(this, x, y);
};

Graph.prototype.dblclick = function(x, y) {
    this.handler.dblclick(this, x, y);
};

Graph.prototype.getPartition = function(callback) {
    var graph = this;
    $.post($SCRIPT_ROOT + "/partition", 
            JSON.stringify({"poly": this.poly, "constrs": this.constrs}),
            function(data){
                graph.partition = data['partition'];
                graph.drawPartitionNames();
                if (callback) {
                    callback();
                }
            }, "json");
};

Graph.prototype.saveSelected = function() {
    if (this.selected != undefined) {
        var pold = this.partition[this.selected];
        pold.A = $("#Coefficients_Matrix").val();
        pold.b_space = this.par_canvas.poly;
    }
};

Graph.prototype.select = function(i) {
    this.saveSelected();
    this.selected = i;
    var p = this.partition[i];
    $("#System_Mode").val(p.name);
    $("#Coefficients_Matrix").val(p.A);
    this.par_canvas.erase();
    this.par_canvas.poly = p.b_space;
    this.par_canvas.drawPoly(p.b_space);
    if (p.b_space_synth) {
        for (var i = 0; i < p.b_space_synth.length; i++) {
            this.par_canvas.drawPoly(p.b_space_synth[i], true);
        }
    }
};

Graph.prototype.selectGraph = function(i) {
    $("#graph" + this.selected_graph).attr("hidden", "true");
    $("#graph" + i).removeAttr("hidden");
    this.selected_graph = i;
};

Graph.prototype.synthesize = function() {
    var m = this.getModel();
    var graph = this;

    $.post($SCRIPT_ROOT + "/synthesize", 
            JSON.stringify(m),
            function(data){
                res = data['result'];
                psets = res['psets']
                for (var i = 0; i < graph.partition.length; i++) {
                    var p = graph.partition[i];
                    p.b_space_synth = psets[i];
                }
                for (var i = 0; i < res['imgs'].length; i++) {
                    img = "<img src={0}/get_image?file={1} id=graph{2} hidden=true/>"
                    $("#graphs").append(img.format(
                                $SCRIPT_ROOT, res['imgs'][i], i));
                }
                graph.selected_graph = 0;
                graph.selectGraph(0);
                console.log(data)
            }, "json");

};

Graph.prototype.getModel = function() {
    var init = $("#Initial_State").val();
    var spec = $("#Spec").val();
    if (this.partition != undefined) {
        this.saveSelected();
        var ps = [];
        for (var i = 0; i < this.partition.length; i++) {
            var p = this.partition[i];
            ps.push({
                A: p.A,
                b_space: p.b_space
            });
        }
        var ts = {
            constrs: this.constrs,
            poly: this.poly,
            pars: ps,
            init: init,
            spec: spec
        };
        return ts;
    }
};

Graph.prototype.loadModel = function(ts) {
    this.poly = ts.poly;
    this.poly_closed = true;
    this.constrs = ts.constrs;
    this.renderPoly();
    this.renderConstrs();
    $("#Initial_State").val(ts.init);
    $("#Spec").val(ts.spec);
    var graph = this;
    this.getPartition(function() {
        for (var i = 0; i < ts.pars.length; i++) {
            var p = graph.partition[i];
            p.A = ts.pars[i].A;
            p.b_space = ts.pars[i].b_space;
        }
    });
}

function addPolygonPoint(graph, x, y) {
    if (graph.poly == undefined) {
        graph.poly = [];
        graph.poly_closed = false;
    }
    if (!graph.poly_closed) {
        graph.poly.push([x,y]);
        graph.renderPoly();
    }
}

function closePolygon(graph, x, y) {
    if (!graph.poly_closed && graph.poly.length >= 2) {
        graph.poly.push([x, y]);
        graph.poly_closed = true;
        graph.renderPoly();
    }
}

function addConstrPoint(graph, x, y) {
    if (graph.cur_constr.length == 0) {
        graph.cur_constr.push([x, y]);
    } else {
        graph.cur_constr.push([x, y]);
        graph.renderCurConstr();
        graph.constrs.push(graph.cur_constr);
        graph.cur_constr = []
    }

}

function selectSet(graph, x, y) {
    parts = graph.partition;
    for (var i = 0; i < parts.length; i++) {
        if (contains(parts[i].constrs, x, y)) {
            graph.select(i);
        }
    }
}

function contains(constrs, x, y) {
    for (var i = 0; i < constrs.length; i++) {
        c = constrs[i];
        if (c[0] + x * c[1] + y * c[2] <= 0) {
            return false;
        }
    }
    return true;
}

function Handler(click, dblclick) {
    this.click = click;
    this.dblclick = dblclick;
}

var Handlers = {
    "poly": new Handler(addPolygonPoint, closePolygon),
    "const": new Handler(addConstrPoint, function(graph, x, y){}),
    "select": new Handler(selectSet, function(graph, x, y){})
}

function getPos(e) {
    var x;
    var y;
    if (e.pageX || e.pageY) { 
    x = e.pageX;
    y = e.pageY;
    }
    else { 
    x = e.clientX + document.body.scrollLeft + document.documentElement.scrollLeft; 
    y = e.clientY + document.body.scrollTop + document.documentElement.scrollTop; 
    } 
    x -= e.target.offsetLeft;
    y -= e.target.offsetTop;
    return [x, y];
}

function handleClick(event) {
    var p = getPos(event);
    var x = p[0], y = p[1];

    var graph = event.data;
    var p = graph.toScale(x, y);
    graph.click(p[0], p[1]);
}

function handleDblClick(event) {
    var p = getPos(event);
    var x = p[0], y = p[1];

    var graph = event.data;
    var p = graph.toScale(x, y);
    graph.dblclick(p[0], p[1]);
}

function draw() {
    var myGraph = new Graph({
        canvasId: 'canvas',
        minX: -10,
        minY: -10,
        maxX: 10,
        maxY: 10,
        unitsPerTick: 1,
        main: true
    });
    myGraph.handler = Handlers["poly"];
    $("input[name='mode']").click(function(event) {
        myGraph.handler = Handlers[event.target.value];
    });
    $("#partition").click(function(event) {
        myGraph.getPartition();
        $("input[name=mode][value=select]").click();
    });
    $("#synthesize").click(function(event) {
        myGraph.synthesize();
    });
    $("#save").click(function(event) {
        $("#model").val(JSON.stringify(myGraph.getModel()));
    });
    $("#load").click(function(event) {
        myGraph.loadModel(JSON.parse($("#model").val()));
    });
    $("#model").load($EXAMPLE_MODEL);

    window.myGraph = myGraph;
}


$(draw);

