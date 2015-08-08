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

    // constants
    this.axisColor = '#aaa';
    this.font = '8pt Calibri';
    this.tickSize = 20;

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
    this.drawXAxis();
    this.drawYAxis();
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
        context.fillText(unit, xPos, this.centerY + this.tickSize / 2 + 3);
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
        context.fillText(unit, xPos, this.centerY + this.tickSize / 2 + 3);
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
        context.fillText(unit, this.centerX - this.tickSize / 2 - 3, yPos);
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
        context.fillText(unit, this.centerX - this.tickSize / 2 - 3, yPos);
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

Graph.prototype.transformPoint = function(x, y) {
    return [(x - this.centerX) / this.scaleX, - (y - this.centerY) / this.scaleY]
};

Graph.prototype.renderPoly = function() {
    var context = this.context;
    context.save();
    this.transformContext();
    context.lineWidth = 0.1;

    if (this.poly.length > 1) {
        context.beginPath();
        context.moveTo(this.poly[0][0], this.poly[0][1]);
        for (var i = 1; i < this.poly.length; i++) {
            context.lineTo(this.poly[i][0], this.poly[i][1]);
        }
        if (this.poly_closed) {
            context.closePath();
        }
        context.stroke();
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

}

Graph.prototype.click = function(x, y) {
    this.handler.click(this, x, y);
};

Graph.prototype.dblclick = function(x, y) {
    this.handler.dblclick(this, x, y);
};

function addPolygonPoint(graph, x, y) {
    if (!graph.poly_closed) {
        graph.poly.push([x,y]);
        graph.renderPoly();
    }
}

function closePolygon(graph, x, y) {
    if (!graph.poly_closed) {
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

function Handler(click, dblclick) {
    this.click = click;
    this.dblclick = dblclick;
}

var Handlers = {
    "poly": new Handler(addPolygonPoint, closePolygon),
    "const": new Handler(addConstrPoint, function(x, y){})
}

function handleClick(event) {
    var x = event.clientX - event.target.offsetLeft;
    var y = event.clientY - event.target.offsetTop;

    var graph = event.data;
    var p = graph.transformPoint(x, y);
    graph.click(p[0], p[1])
}

function handleDblClick(event) {
    var x = event.clientX - event.target.offsetLeft;
    var y = event.clientY - event.target.offsetTop;

    var graph = event.data;
    var p = graph.transformPoint(x, y);
    graph.dblclick(p[0], p[1])
}

function draw() {
    var myGraph = new Graph({
        canvasId: 'canvas',
        minX: -10,
        minY: -10,
        maxX: 10,
        maxY: 10,
        unitsPerTick: 1
    });
    myGraph.handler = Handlers["poly"];
    $("input[name='mode']").click(function(event) {
        myGraph.handler = Handlers[event.target.value];
    });
}


$(draw);

