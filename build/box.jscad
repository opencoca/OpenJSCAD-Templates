
function main(params) {
  var config = params.config;
  var plates = {};

  // set render quality of dxf export
  CSG.defaultResolution2D = 32;
  
  if (config.render.shell) {
    plates = Object.assign(plates, Shell.getPlates(config));
  }
  
  plates = Object.assign(plates, Shelves.getPlates(config));
 
  // shell finger joints for itself
  if (plates.shell) {
    plates.shell = Joints.fingerJoint(plates.shell, Utils.getFingerJointOptions(config));
  }

  // add holes to side shells to mount horizontal shelf dividers
  // TODO: optimize
  if (plates.shell && plates.shelves) {
    for (var i=0; i < plates.shell.length; i++) {
      for (var j=0; j < plates.shelves.length; j++) {
       plates.shell[i] = difference(plates.shell[i], plates.shelves[j]);
      }
      for (var j=0; j < plates.dividers.length; j++) {
       plates.shell[i] = difference(plates.shell[i], plates.dividers[j]);
      }
    }
  }

  return renderAllPlates(config, plates);
}

function renderAllPlates(config, plates) {
  var plateArray = [];

  for (var plateKey in plates) {
    plateArray = plateArray.concat(plates[plateKey])
  }

  if (config.buildType === 'DXF') {
    var plates2d = []; 
    var offset = 0;

    for (var i=0; i<plateArray.length; i++) {
      // scale all plates down by 10
      //plateArray[i] = plateArray[i].scale(0.1);

      if (i > 0) {
        var bounds1 = Plate.CSGToCAG(plateArray[i-1]).center().getBounds();
        var bounds2 = Plate.CSGToCAG(plateArray[i]).center().getBounds();
        offset += ((bounds1[1].x - bounds1[0].x) / 2) + ((bounds2[1].x - bounds2[0].x) / 2) + 10
      }
      
      plates2d.push(Plate.CSGToCAG(plateArray[i]).center().translate([offset, 0]));
    }

    return plates2d;
  }
  else {
    return plateArray;
  }
}

function getParameterDefinitions()
{
  var params = [];
  return params;
}



// author     : Joost Nieuwenhuijse
// license    : MIT License

var Plate = {
  CSGToCAG: function(plate) {
    if(!("platebasis" in plate.properties))
    {
      throw new Error("Plates should be created using getStockPlate()");
    }
    var plate2d = plate.projectToOrthoNormalBasis(plate.properties.platebasis);
    return plate2d;
  },

  CAGToCSG: function(plate2d, platebasis, thickness) {
    var basisinversematrix = platebasis.getInverseProjectionMatrix();
    var plate_reprojected = plate2d.extrude({offset: [0,0,thickness]}).translate([0,0,-thickness/2]);
    plate_reprojected = plate_reprojected.transform(basisinversematrix);
    plate_reprojected.properties.platebasis = platebasis;
    return plate_reprojected;
  },

  fix: function(config, plate) {
    return Plate.CAGToCSG(Plate.CSGToCAG(plate), plate.properties.platebasis, config.thickness); 
  },

  removeWithNormal: function(plates, normalvector) {
    normalvector = new CSG.Vector3D(normalvector);
    var result = [];
    plates.map(function(plate){
      if(!("platebasis" in plate.properties))
      {
        throw new Error("Plates should be created using getStockPlate()");
      } 
      if(plate.properties.platebasis.plane.normal.dot(normalvector) < 0.9999)
      {
        result.push(plate);
      } 
    });
    return result;
  },

  getStock: function(width, height, thickness) {
    var result = CSG.cube({radius: [width/2, height/2, thickness/2]});
    result.properties.platebasis = CSG.OrthoNormalBasis.Z0Plane();
    return result.translate([width/2, height/2, thickness/2]);
  }
};
var _cache = {};

var Utils = {
  // get total height of all shelves
  getShelvesHeight: function(config) {
    if (_cache.getShelvesHeight) return _cache.getShelvesHeight;
    var height = 0;

    for (var i=0; i < config.shelfCount; i++) {
      // add shelf horizontal divider thicknesses
      if (i > 0) {
        height += config.thickness;
      }

      height += config.shelves['Shelf ' + (i+1)].height;
    }

    _cache.getShelvesHeight = height;

    return height;
  },

  getFitShelfCount(config) {
    var count = 0;

    for (var i=0; i < config.shelfCount; i++) {
      var shelf = config.shelves['Shelf ' + (i + 1)];

      if (shelf.height === 0) {
        count++;
      }
    }

    return count;
  },

  getAvailableFitShelfSpaceHeight(config) {
    var height = Utils.shelfAreaHeight(config);

    // subtract horizontal divider heights from available space
    height -= (config.shelfCount-1) * config.thickness;

    for (var i=0; i < config.shelfCount; i++) {
      var shelf = config.shelves['Shelf ' + (i + 1)];

      // non auto shelf height, subtract from available fit space
      if (shelf.height > 0) {
        height -= shelf.height;
      }
    }

    return height;
  },

  getShelfInnerHeight: function(config, shelfIndex) {
    var shelfConfig = config.shelves['Shelf ' + (shelfIndex+1)];

    // fit
    if (shelfConfig.height === 0) {
      return Utils.getAvailableFitShelfSpaceHeight(config) / Utils.getFitShelfCount(config);
    }

    return shelfConfig.height;
  },

  getAvailableFitModuleSpace(config, shelfConfig) {
    var width = Utils.shelfAreaWidth(config);
    width -= (shelfConfig.moduleCount-1) * config.thickness;

    for (var i=0; i < shelfConfig.moduleCount; i++) {
      var module = shelfConfig.modules['Module ' + (i + 1)];

      if (module.width > 0) {
        width -= module.width;
      }
    }

    return width;
  },

  getFitModuleCount(config, shelfConfig) {
    var count = 0;

    for (var i=0; i < shelfConfig.moduleCount; i++) {
      var module = shelfConfig.modules['Module ' + (i + 1)];

      if (module.width === 0) {
        count++;
      }
    }

    return count;
  },

  getModuleInnerWidth: function(config, shelfConfig, moduleIndex) {
    var moduleConfig = shelfConfig.modules['Module ' + (moduleIndex+1)];

    // fit
    if (moduleConfig.width === 0) {
      return Utils.getAvailableFitModuleSpace(config, shelfConfig) / Utils.getFitModuleCount(config, shelfConfig);
    }

    return moduleConfig.width;
  },

  // https://stackoverflow.com/questions/6211613/testing-whether-a-value-is-odd-or-even
  isEven: function(n) {
    return n == parseFloat(n)? !(n%2) : void 0;
  },

  getFingerJointOptions: function(config) {
    return {
      margin: 0,
      cutterRadius: 0,
      fingerWidth: config.fingerJointWidth
    };
  },

  boxHeight: function(config) {
    return config.boxHeight;
  },

  boxWidth: function(config) {
    return config.boxWidth;
  },

  boxDepth: function(config) {
    return config.boxDepth;
  },

  shelfAreaHeight: function(config) {
    return Utils.boxHeight(config) - (config.thickness*2);
  },

  shelfAreaWidth: function(config) {
    return Utils.boxWidth(config) - (config.thickness*2);
  },

  shelfAreaDepth: function(config) {
    return Utils.boxDepth(config) - (config.thickness);
  }
};


var Joints = {
  fingerJointAdd: function(plates, newplate, options) {
    var result = plates.slice(0);
    var numplates = plates.length;
    for(var plateindex1 = 0; plateindex1 < numplates; plateindex1++)
    {
      var joined = Joints.fingerJointTwo(result[plateindex1], newplate, options);
      result[plateindex1] = joined[0];
      newplate = joined[1];
    }
    result.push(newplate);
    return result;
  },

  addFingers: function(config, plates, options) {
    var fingers = [];

    var plate = plates[0];
    var bounds = plate.getBounds();
    
    if (options.where === "TOP") {
      var length = (bounds[1].x - bounds[0].x);
      var jointCount = ((Math.ceil((length-(config.thickness*4)) / config.fingerJointWidth)/2)*2);
      var offsetX = (length - (jointCount * config.fingerJointWidth))/2;

      // center odd finger joint coints
      if (Utils.isEven(jointCount)) {
        offsetX += config.fingerJointWidth/2;
      }

      for (var j=0; j < jointCount; j += 2) {
        var finger = cube([config.fingerJointWidth, config.thickness, config.thickness])
          .translate([j*config.fingerJointWidth, -config.thickness, 0])
          .translate([offsetX, 0, 0]);

        plate = union(plate, finger);
      }
    }
    else if (options.where === "RIGHT") {
      var length = (bounds[1].y - bounds[0].y);
      var jointCount = ((Math.floor((length-(config.thickness*4)) / config.fingerJointWidth)/2)*2);
      var offsetY = (length - (jointCount * config.fingerJointWidth))/2;

      // center odd finger joint coints
      if (!Utils.isEven(jointCount)) {
        offsetY -= config.fingerJointWidth/2;
      }

      for (var j=0; j < jointCount; j += 2) {
        if (options.skip_even && Utils.isEven(j/2)) continue;
        if (options.skip_odd && !Utils.isEven(j/2)) continue;
        var finger = cube([config.thickness, config.fingerJointWidth, config.thickness])
          .translate([-config.thickness, j*config.fingerJointWidth, 0])
          .translate([0, offsetY, 0]);

        plate = union(plate, finger);
      }
    }
    else if (options.where === "LEFT") {
      var length = (bounds[1].y - bounds[0].y);
      var jointCount = ((Math.floor((length-(config.thickness*4)) / config.fingerJointWidth)/2)*2);
      var offsetY = (length - (jointCount * config.fingerJointWidth))/2;

      // center odd finger joint coints
      if (!Utils.isEven(jointCount)) {
        offsetY -= config.fingerJointWidth/2;
      }

      for (var j=0; j < jointCount; j += 2) {
        if (options.skip_even && Utils.isEven(j/2)) continue;
        if (options.skip_odd && !Utils.isEven(j/2)) continue;
        var finger = cube([config.thickness, config.fingerJointWidth, config.thickness])
          .translate([-config.thickness, j*config.fingerJointWidth, 0])
          .translate([(bounds[1].x - bounds[0].x), offsetY, 0]);

        plate = union(plate, finger);
      }
    }

    return Plate.fix(config, plate);
  },

  // Finger joint between multiple plates:
  fingerJoint: function(plates, options) {
    var result = plates.slice(0);
    var numplates = plates.length;
    var maxdelta = Math.floor(numplates/2);
    for(var delta=1; delta <= maxdelta; delta++)
    { 
      for(var plateindex1 = 0; plateindex1 < numplates; plateindex1++)
      {
        var plateindex2 = plateindex1 + delta;
        if(plateindex2 >= numplates) plateindex2 -= numplates; 
        
        var joined = Joints.fingerJointTwo(result[plateindex1], result[plateindex2], options);
        result[plateindex1] = joined[0];
        result[plateindex2] = joined[1];
        if(delta*2 >= numplates)
        {
          // numplates is even
          if(plateindex1*2 >= numplates)
          {
            // and we've done the first half: we're done
            break;
          }
        }
      }
    }
    return result;
  },

  fingerJointTwo: function(plate1, plate2, options) {
    if(!options) options = {};
    if(!("platebasis" in plate1.properties))
    {
      throw new Error("Plates should be created using getStockPlate()");
    } 
    if(!("platebasis" in plate2.properties))
    {
      throw new Error("Plates should be created using getStockPlate()");
    } 
    // get the intersection solid of the 2 plates:
    var intersection = plate1.intersect(plate2);
    if(intersection.polygons.length === 0)
    {
      // plates do not intersect. Return unmodified:
      return [plate1, plate2];
    }
    else
    {
      var plane1 = plate1.properties.platebasis.plane; 
      var plane2 = plate2.properties.platebasis.plane;
      // get the intersection line of the 2 center planes:
      var jointline = plane1.intersectWithPlane(plane2);
      // Now we need to find the two endpoints on jointline (the points at the edges of intersection):
      // construct a plane perpendicular to jointline:
      plane1 = CSG.Plane.fromNormalAndPoint(jointline.direction, jointline.point);
      // make the plane into an orthonormal basis:
      var basis1 = new CSG.OrthoNormalBasis(plane1);
      // get the projection matrix for the orthobasis:
      var matrix = basis1.getProjectionMatrix();
      // now transform the intersection solid:
      var intersection_transformed = intersection.transform(matrix);
      var bounds = intersection_transformed.getBounds();
      // now we know the two edge points. The joint line runs from jointline_origin, in the 
      // direction jointline_direction and has a length jointline_length (jointline_length >= 0) 
      var jointline_origin = jointline.point.plus(jointline.direction.times(bounds[0].z));
      var jointline_direction = jointline.direction;
      var jointline_length = (bounds[1].z - bounds[0].z); 
    
      var fingerwidth = options.fingerWidth || (jointline_length / 4); 
      var numfingers=Math.round(jointline_length / fingerwidth);
      if(numfingers < 2) numfingers=2;
      fingerwidth = jointline_length / numfingers;
      
      var margin = options.margin || 0;
      var cutterRadius = options.cutterRadius || 0; 
      var results = [];
      for(var plateindex = 0; plateindex < 2; plateindex++)
      {
        var thisplate = (plateindex == 1)? plate2:plate1;
        var otherplate = (plateindex == 1)? plate1:plate2;
        // create a new orthonormal basis for this plate, such that the joint line runs in the positive x direction:
        var platebasis = new CSG.OrthoNormalBasis(thisplate.properties.platebasis.plane, jointline_direction);
        // get the 2d shape of our plate:
        var plate2d = thisplate.projectToOrthoNormalBasis(platebasis);
        var jointline_origin_2d = platebasis.to2D(jointline_origin);
        matrix = platebasis.getProjectionMatrix();
        intersection_transformed = intersection.transform(matrix);
        bounds = intersection_transformed.getBounds();
        var maxz = bounds[1].z;
        var minz = bounds[0].z;
        var maxy = bounds[1].y + margin/2; 
        var miny = bounds[0].y - margin/2;
        
        var cutouts2d = [];
        var othercutouts2d = [];
        for(var fingerindex = 0; fingerindex < numfingers; fingerindex++)
        {
          if( (plateindex === 0) && ((fingerindex & 1)===0) ) continue;
          if( (plateindex  == 1) && ((fingerindex & 1)!==0) ) continue;
          var minx = jointline_origin_2d.x + fingerindex * fingerwidth - margin/2;
          var maxx = minx + fingerwidth + margin;
          var cutout = Joints.createRectCutoutWithCutterRadius(minx, miny, maxx, maxy, cutterRadius, plate2d);
          //fingerindex == 0 || fingerindex == numfingers-1 ? othercutouts2d.push(cutout) : cutouts2d.push(cutout)
          cutouts2d.push(cutout);
        }
        var cutout2d = new CAG().union(cutouts2d);
        var cutout3d = cutout2d.extrude({offset: [0,0,maxz-minz]}).translate([0,0,minz]);
        cutout3d = cutout3d.transform(platebasis.getInverseProjectionMatrix());
        var thisplate_modified = thisplate.subtract(cutout3d);
        results[plateindex] = thisplate_modified;  
      } 
      return results;
    }
  },

  // Create a rectangular cutout in 2D
  // minx, miny, maxx, maxy: boundaries of the rectangle
  // cutterRadius: if > 0, add extra cutting margin at the corners of the rectangle
  // plate2d is the 2d shape from which the cutout will be subtracted
  // it is tested at the corners of the cutout rectangle, to see if do need to add the extra margin at that corner
  createRectCutoutWithCutterRadius: function(minx, miny, maxx, maxy, cutterRadius, plate2d) {
    var deltax = maxx-minx;
    var deltay = maxy-miny;  
    var cutout = CAG.rectangle({radius: [(maxx-minx)/2, (maxy-miny)/2], center: [(maxx+minx)/2, (maxy+miny)/2]});
    var cornercutouts = [];
    if(cutterRadius > 0)
    {
      var extracutout = cutterRadius * 0.2;
      var hypcutterradius = cutterRadius / Math.sqrt(2.0);
      var halfcutterradius = 0.5 * cutterRadius;
      var dcx, dcy;
      if(deltax > 3*deltay)
      {
        dcx = cutterRadius + extracutout/2;
        dcy = extracutout / 2;
      }
      else if(deltay > 3*deltax)
      {
        dcx = extracutout / 2;
        dcy = cutterRadius + extracutout/2;
      }
      else
      {
        dcx = hypcutterradius-extracutout/2;
        dcy = hypcutterradius-extracutout/2;
      }
      for(var corner = 0; corner < 4; corner++)
      {
        var cutoutcenterx = (corner & 2)? (maxx-dcx):(minx+dcx);
        var cutoutcentery = (corner & 1)? (maxy-dcy):(miny+dcy);
        var cornercutout = CAG.rectangle({radius: [cutterRadius+extracutout/2, cutterRadius+extracutout/2], center: [cutoutcenterx, cutoutcentery]});
        var testrectacenterx = (corner & 2)? (maxx-halfcutterradius):(minx+halfcutterradius);
        var testrectbcenterx = (corner & 2)? (maxx+halfcutterradius):(minx-halfcutterradius);
        var testrectacentery = (corner & 1)? (maxy+halfcutterradius):(miny-halfcutterradius);
        var testrectbcentery = (corner & 1)? (maxy-halfcutterradius):(miny+halfcutterradius);
        var testrecta = CAG.rectangle({radius: [halfcutterradius, halfcutterradius], center: [testrectacenterx, testrectacentery]}); 
        var testrectb = CAG.rectangle({radius: [halfcutterradius, halfcutterradius], center: [testrectbcenterx, testrectbcentery]});
        if( (plate2d.intersect(testrecta).sides.length > 0)  &&
         (plate2d.intersect(testrectb).sides.length > 0) )
        {
          cornercutouts.push(cornercutout);
        }
      }
    }
    if(cornercutouts.length > 0)
    {
      cutout = cutout.union(cornercutouts);
    } 
    return cutout;
  }
};
var Shell = {
  getPlates: function(config) {
    var boxWidth = Utils.boxWidth(config);
    var boxHeight = Utils.boxHeight(config);
    var boxDepth = Utils.boxDepth(config);

    var plates = [
      Shell.getTopPlate(config, boxWidth, boxHeight, boxDepth),
      Shell.getLeftSidePlate(config, boxWidth, boxHeight, boxDepth),
      Shell.getRightSidePlate(config, boxWidth, boxHeight, boxDepth),
      Shell.getBottomPlate(config, boxWidth, boxHeight, boxDepth),
      Shell.getBackPlate(config, boxWidth, boxHeight, boxDepth)
    ];

    return {
      shell: plates
    };
  },

  getTopPlate: function(config, boxWidth, boxHeight, boxDepth) {
    var topPlate = Plate.getStock(boxWidth, boxDepth, config.thickness)
      .translate([0,0,boxHeight-config.thickness]);
    
    return topPlate.setColor([198/255, 144/255, 35/255]);
  },

  getLeftSidePlate: function(config, boxWidth, boxHeight, boxDepth) {
    var sidePlate = Plate.getStock(boxHeight, boxDepth, config.thickness)
      .rotateY(90)
      .translate([boxWidth-config.thickness,0,boxHeight]);

    return sidePlate.setColor([122/255, 89/255, 23/255]);
  },

  getRightSidePlate: function(config, boxWidth, boxHeight, boxDepth) {
    var sidePlate = Plate.getStock(boxHeight, boxDepth, config.thickness)
      .rotateY(90)
      .translate([0,0,boxHeight]);

    return sidePlate.setColor([122/255, 89/255, 23/255]);
  },

  getBottomPlate: function(config, boxWidth, boxHeight, boxDepth) {
    var bottomPlate = Plate.getStock(boxWidth, boxDepth, config.thickness);
    
    return bottomPlate.setColor([198/255, 144/255, 35/255]);
  },

  getBackPlate: function(config, boxWidth, boxHeight, boxDepth) {
    var backPlate = Plate.getStock(boxWidth, boxHeight, config.thickness)
      .rotateX(90)
      .translate([0,config.thickness,0]);

      return backPlate.setColor([170/255, 95/255, 33/255]);
  }
};
var Drawer = {
  getPlates: function(config, shelfIndex, moduleIndex) {
    var _m = config.moduleMargin;
    var shelfConfig = config.shelves['Shelf ' + (shelfIndex+1)];
    var moduleConfig = shelfConfig.modules['Module ' + (moduleIndex+1)];
    var moduleWidth = Utils.getModuleInnerWidth(config, shelfConfig, moduleIndex) - (_m*2);
    var moduleHeight = Utils.getShelfInnerHeight(config, shelfIndex) - (_m*2);
    var moduleDepth = Utils.shelfAreaDepth(config) - (_m);

    var plates = [
      Drawer.getBottomPlate(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth),
      Drawer.getRightPlate(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth),
      Drawer.getLeftPlate(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth),
      Drawer.getBackPlate(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth),
      Drawer.getFrontPlate(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth)
    ];

    plates = Joints.fingerJoint(plates, Utils.getFingerJointOptions(config));

    return {
      drawer: plates
    };
  },

  getPullTab(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth, tabDepth) {
    var tab = cube({size: [moduleWidth, tabDepth-config.thickness, config.thickness]})
    var r = config.thickness;

    // dem rounded corners tho
    tab = union(
      tab,
      cylinder({r: r, h: config.thickness}).translate([r,tabDepth-r,0]),
      cylinder({r: r, h: config.thickness}).translate([moduleWidth - r, tabDepth-r, 0]),
      cube({size: [moduleWidth - (r*2), r, config.thickness]}).translate([r, tabDepth-r,config.thickness])
    );

    return tab;
  },

  getBottomPlate: function(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth) {
    var tabDepth = config.thickness*2;
    var plate = Plate.getStock(moduleWidth, moduleDepth, config.thickness)
      .setColor([111/255,0,0]);

    var tab = Drawer.getPullTab(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth, tabDepth)
      .translate([0,moduleDepth,0])

    plate = union(plate, tab)

    plate = Plate.fix(config, plate);

    return plate;
  },

  getRightPlate: function(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth) {
    var plate = Plate.getStock(moduleHeight, moduleDepth, config.thickness)
      .rotateY(90)
      .translate([0,0,moduleHeight])
      .setColor([155/255,0,0]);

    return plate;
  },

  getLeftPlate: function(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth) {
    var plate = Plate.getStock(moduleHeight, moduleDepth, config.thickness)
      .rotateY(90)
      .translate([moduleWidth-config.thickness,0,moduleHeight])
      .setColor([155/255,0,0]);

    return plate;
  },

  getBackPlate: function(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth) {
    var plate = Plate.getStock(moduleWidth, moduleHeight, config.thickness);
    var r = Math.min(moduleHeight, moduleWidth) / 3;

    plate = difference(plate, 
      cylinder({r: r, h: config.thickness*2}).translate([moduleWidth/2,moduleHeight + (r/2.5),-config.thickness]));

    plate = plate.rotateX(90)
      .translate([0,config.thickness,0])
      .setColor([200/255,0,0]);

    return plate;
  },

  getFrontPlate: function(config, moduleConfig, moduleWidth, moduleHeight, moduleDepth) {
    var plate = Plate.getStock(moduleWidth, moduleHeight, config.thickness);

    plate = plate.rotateX(90)
      .translate([0,moduleDepth,0])
      .setColor([200/255,0,0]);

    return plate;
  }
};
var Shelves = {
  getPlates: function(config) {
    var shelfAreaDepth = Utils.shelfAreaDepth(config);
    var shelfAreaWidth = Utils.shelfAreaWidth(config);
    var shelfAreaHeight = Utils.shelfAreaHeight(config);

    var plates = {
      shelves: Shelves.renderHorizontal(config, shelfAreaWidth, shelfAreaDepth, shelfAreaHeight),
      dividers: Shelves.renderVertical(config, shelfAreaWidth, shelfAreaDepth, shelfAreaHeight),
      modules: Shelves.renderModules(config, shelfAreaWidth, shelfAreaDepth, shelfAreaHeight)
    }

    if (plates.dividers && plates.shelves) {
      for (var i=0; i < plates.shelves.length; i++) {
        for (var j=0; j < plates.dividers.length; j++) {
          plates.shelves[i] = difference(plates.shelves[i], plates.dividers[j]);
        }
      }
    }

    if (!config.render.shelves) {
      delete plates.shelves;
      delete plates.dividers;
    }
    if (!config.render.modules) {
      delete plates.modules;
    }

    return plates;
  },

  renderHorizontal(config, shelfAreaWidth, shelfAreaDepth, shelfAreaHeight) {
    var plates = [];
    var offsetZ = shelfAreaHeight + (config.thickness);

    for (var i=0; i < config.shelfCount; i++) {
      var horizontalDivider = Plate.getStock(shelfAreaWidth, shelfAreaDepth, config.thickness);
    
      horizontalDivider = horizontalDivider;

      if (i > 0) {
        horizontalDivider = Joints.addFingers(config, [horizontalDivider], {where:'TOP'})
        horizontalDivider = Joints.addFingers(config, [horizontalDivider], {where:'RIGHT'})
        horizontalDivider = Joints.addFingers(config, [horizontalDivider], {where:'LEFT'})

        horizontalDivider = horizontalDivider.setColor([96/255, 117/255, 150/255]);

        plates.push(horizontalDivider.translate([config.thickness, config.thickness, offsetZ]));
      }

      offsetZ -= Utils.getShelfInnerHeight(config, i);
      offsetZ -= config.thickness;
    }

    return plates;
  },

  renderVertical(config, shelfAreaWidth, shelfAreaDepth, shelfAreaHeight) {
    var plates = [];
    var offsetZ = shelfAreaHeight;

    for (var i=0; i < config.shelfCount; i++) {
      var shelfConfig = config.shelves['Shelf ' + (i+1)];
      var offsetX = shelfAreaWidth;

      for (var j=0; j < shelfConfig.moduleCount; j++) {
        var moduleConfig = shelfConfig.modules['Module ' + (j+1)];
        var shelfInnerHeight = Utils.getShelfInnerHeight(config, i);
        var verticalDivider = Plate.getStock(shelfInnerHeight, shelfAreaDepth, config.thickness);

        verticalDivider = Joints.addFingers(config, [verticalDivider], {where:'TOP'})
        verticalDivider = Joints.addFingers(config, [verticalDivider], {where:'RIGHT', skip_odd: true})
        verticalDivider = Joints.addFingers(config, [verticalDivider], {where:'LEFT', skip_even: true})
        
        verticalDivider = verticalDivider.rotateY(90)
          .translate([offsetX, config.thickness, offsetZ+config.thickness]);
      
        verticalDivider = verticalDivider.setColor([157/255, 189/255, 242/255]);

        offsetX -= Utils.getModuleInnerWidth(config, shelfConfig, j);

        if (j > 0) {
          offsetX -= config.thickness;
          plates.push(verticalDivider)
        }
      }

      offsetZ -= Utils.getShelfInnerHeight(config, i);
      offsetZ -= config.thickness;
    }

    return plates;
  },

  renderModules(config, shelfAreaWidth, shelfAreaDepth, shelfAreaHeight) {
    var plates = [];
    var offsetZ = shelfAreaHeight + config.thickness;
    var _m = config.moduleMargin;

    for (var i=0; i < config.shelfCount; i++) {
      var shelfConfig = config.shelves['Shelf ' + (i+1)];
      var offsetX = shelfAreaWidth + config.thickness;

      offsetZ -= Utils.getShelfInnerHeight(config, i);

      for (var j=0; j < shelfConfig.moduleCount; j++) {
        var moduleConfig = shelfConfig.modules['Module ' + (j+1)];

        offsetX -= Utils.getModuleInnerWidth(config, shelfConfig, j);

        if (moduleConfig.type === "Drawer") {
          var drawer = Drawer.getPlates(config, i, j);

          for (var k=0; k < drawer.drawer.length; k++) {
            drawer.drawer[k] = drawer.drawer[k]
              .translate([offsetX + (_m), config.thickness + (_m), offsetZ + (_m)])
          }

          plates = plates.concat(drawer.drawer);
        }

        offsetX -= config.thickness;
      }

      offsetZ -= config.thickness;
    }

    return plates;
  }
}
