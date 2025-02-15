
 
 /*
COMP 3490 Assignment 1 Template - Math Library
*/

// Helper function to check if a point is inside a triangle using barycentric coordinates
boolean pointInTriangle(float x, float y, PVector[] v) {
  float[] bary = computeBarycentricCoordinates(x, y, v);
  float u = bary[0], vCoord = bary[1], w = bary[2];
  // Allow a smallish epsilon to check 4 edge cases due to floating point precision
  float epsilon = -0.01f;
  return (u >= epsilon) && (vCoord >= epsilon) && (w >= epsilon);
}

// Computate barycentric coordinates for a point i.e x, y with respect to triangle vertices
float[] computeBarycentricCoordinates(float x, float y, PVector[] v) {
  float x0 = v[0].x, y0 = v[0].y;
  float x1 = v[1].x, y1 = v[1].y;
  float x2 = v[2].x, y2 = v[2].y;

  float denom = (y1 - y2)*(x0 - x2) + (x2 - x1)*(y0 - y2);
  if (abs(denom) < 1e-6) {
    // this is Degenerate triangle
    return new float[] { -1, -1, -1 };
  }

  float u = ((y1 - y2)*(x - x2) + (x2 - x1)*(y - y2)) / denom;
  float vCoord = ((y2 - y0)*(x - x2) + (x0 - x2)*(y - y2)) / denom;
  float w = 1 - u - vCoord;
  return new float[] { u, vCoord, w };
}

// Function to check if a triangle is front-facing
boolean isFrontFacing(PVector[] projectedVertices) {
  PVector edge1 = PVector.sub(projectedVertices[1], projectedVertices[0]);
  PVector edge2 = PVector.sub(projectedVertices[2], projectedVertices[0]);

  // Cross product in 2D: if z-component is positive, the triangle is front-facing
  float cross = edge1.x * edge2.y - edge1.y * edge2.x;
  return cross >= 0.1; // Allow some small error margin to avoid too much culling
}
