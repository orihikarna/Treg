// https://en.wikibooks.org/wiki/OpenSCAD_User_Manual/Commented_Example_Projects#Icosphere

// Code via reddit with triangle winding fixes, cannot add link due to
// wikibooks considering it spam.

// 4 is the realistic max.
// Don't do 5 or more, takes forever.
// set recursion to the desired level. 0=20 tris, 1=80 tris, 2=320 tris
module icosphere(radius = 10, recursion = 2, icoPnts, icoTris) {
  //t = (1 + sqrt(5))/2;
  //comment from monfera to get verts to unit sphere
  t = sqrt((5 + sqrt(5)) / 10);
  s = sqrt((5 - sqrt(5)) / 10);

  init = (icoPnts || icoTris) ? false : true;//initial call if icoPnts is empty

  // 1 --> draw icosphere from base mesh
  // 2 --> loop through base mesh and subdivide by 4 --> 20 steps
  // 3 --> loop through subdivided mesh and subdivide again (or subdivide by 16) --> 80 steps
  // 4 ...

  verts = [
    [-s, t, 0],//0
    [s, t, 0], 
    [-s, -t, 0], 
    [s, -t, 0], 
    [0, -s, t], 
    [0, s, t], 
    [0, -s, -t], 
    [0, s, -t], 
    [t, 0, -s], 
    [t, 0, s], 
    [-t, 0, -s], 
    [-t, 0, s]
  ];//11

  //base mesh with 20 faces
  tris = [
    //5 faces around point 0
    [0, 5, 11],//0
    [0, 1, 5], 
    [0, 7, 1], 
    [0, 10, 7], 
    [0, 11, 10], 
    // 5 adjacent faces
    [1, 9, 5],//5
    [5, 4, 11], 
    [11, 2, 10], 
    [10, 6, 7], 
    [7, 8, 1], 
    //5 faces around point 3
    [3, 4, 9],//10
    [3, 2, 4], 
    [3, 6, 2], 
    [3, 8, 6], 
    [3, 9, 8], 
    //5 adjacent faces
    [4, 5, 9],//15
    [2, 11, 4], 
    [6, 10, 2], 
    [8, 7, 6], 
    [9, 1, 8]
  ];//19

  if (recursion) {
    verts = (init) ? verts : icoPnts;
    tris = (init) ? tris : icoTris;
    newSegments = recurseTris(verts, tris);
    newVerts = newSegments[0];
    newTris = newSegments[1];
    icosphere(radius, recursion - 1, newVerts, newTris);
  } else if (init) {//draw the base icosphere if no recursion and initial call
    scale(radius)
      polyhedron(verts, tris);
  } else {// if not initial call some recursion has to be happened
    scale(radius)
      polyhedron(icoPnts, icoTris);
  }
}

// Adds verts if not already there,
// takes array of vertices and indices of a tri to expand
// returns expanded array of verts and indices of new polygon with 4 faces
// [[verts],[0,(a),(c)],[1,(b),(a)],[2,(c),(b)],[(a),(b),(c)]]
function addTris(verts, tri) = let(a = getMiddlePoint(verts[tri[0]], verts[tri[1]]),//will produce doubles
b = getMiddlePoint(verts[tri[1]], verts[tri[2]]),//these are unique
c = getMiddlePoint(verts[tri[2]], verts[tri[0]]),//these are unique

aIdx = search(verts, a),//point a already exists
l = len(verts))len(aIdx) ? [concat(verts, 
  [a, b, c]
), 
  [
    [tri[0], l, l + 2],//1
    [tri[1], l + 1, l],//2
    [tri[2], l + 2, l + 1],//3
    [l, l + 1, l + 2]
  ]
] ://4

[concat(verts, 
  [b, c]
), 
  [
    [tri[0], aIdx, l + 1],//1
    [tri[1], l, aIdx],//2
    [tri[2], l + 1, l],//3
    [aIdx, l, l + 1]
  ]
];//4

// Recursive function that does one recursion on the whole icosphere (auto recursion steps derived from len(tris)).
function recurseTris(verts, tris, newTris = [], steps = 0, step = 0) = let(stepsCnt = steps ? steps : len(tris) - 1,//if initial call initialize steps
newSegment = addTris(verts = verts, tri = tris[step]), newVerts = newSegment[0],//all old and new Vertices
newerTris = concat(newTris, newSegment[1])//only new Tris
)(stepsCnt == (step)) ? [newVerts, newerTris] : recurseTris(newVerts, tris, newerTris, stepsCnt, step + 1);

// Get point between two verts on unit sphere.
function getMiddlePoint(p1, p2) = fixPosition((p1 + p2) / 2);

// Fix position to be on unit sphere
function fixPosition(p) = let(l = norm(p))[p.x / l, p.y / l, p.z / l];