/*
COMP 3490 Assignment 1 Template - Main Rendering Code
*/

// Constants variable
final color BLACK = color(0);
final int nPhi = 15;   //  latitude divisions
final int nTheta = 20; //  longitude divisions

final float sphereRadius = 100;  // Radius of the sphere

final int NORMAL_LENGTH = 20;
final float[] FACE_NORMAL_COLOR = {0f, 1f, 1f}; // cyan color
final float[] VERTEX_NORMAL_COLOR = {1f, 0.6f, 0.1f}; // orange i think


// For 1 triangle mode
Triangle[] singleTriangle;
Triangle[] rotatedSingleTriangle;

// For drawing and rotating the 3D shape (sphere was chosen)
Triangle[] surfaceTessellation;
Triangle[] rotatedSurfaceTessellation;

// to make the image rotate - don't change these values
float theta = 0.0;  // rotation angle
float dtheta = 0.01; // rotation speed

PGraphics buffer;

// Setup function
void setup() {
  // RGB values range over 0..1 rather than 0..255
  colorMode(RGB, 1.0f);

  buffer = createGraphics(600, 600);

  // Create a single triangle in 3D space
  PVector v1 = new PVector(-200, -100, 50);  // Vertex 1
  PVector v2 = new PVector(100, -100, 50);   // Vertex 2
  PVector v3 = new PVector(0, 100, 90);      // Vertex 3

  PVector n1 = new PVector(0, 0, 1);
  PVector n2 = new PVector(0, 0, 1);
  PVector n3 = new PVector(0, 0, 1);

  singleTriangle = new Triangle[] { new Triangle(new PVector[] { v1, v2, v3 }, new PVector[] { n1, n2, n3 }) };
  rotatedSingleTriangle = copyTriangleList(singleTriangle);

  // Tessellate the sphere and store it in surfaceTessellation
  surfaceTessellation = tessellateSphere(sphereRadius, nPhi, nTheta);
  rotatedSurfaceTessellation = copyTriangleList(surfaceTessellation);

  printSettings();
}

void settings() {
  size(600, 600); // hard-coded canvas size, same as the buffer
}

/*
You should read this function carefully and understand how it works,
 but you should not need to edit it
 */
void draw() {
  buffer.beginDraw();
  buffer.colorMode(RGB, 1.0f);
  buffer.background(BLACK);

  /*
  CAUTION: none of your functions should call loadPixels() or updatePixels().
   This is already done in the template. Extra calls will probably break things.
   */
  buffer.loadPixels();

  if (doRotate) {
    theta += dtheta;
    if (theta > TWO_PI) {
      theta -= TWO_PI;
    }
  }

  //do not change these blocks: rotation is already set up for you
  if (displayMode == DisplayMode.TEST_LINES) {
    testBresenham();
  } else if (displayMode == DisplayMode.SINGLE_TRIANGLE) {
    rotateTriangles(singleTriangle, rotatedSingleTriangle, theta);
    drawTriangles(rotatedSingleTriangle);
  } else if (displayMode == DisplayMode.SURFACE) {
    rotateTriangles(surfaceTessellation, rotatedSurfaceTessellation, theta);
    drawTriangles(rotatedSurfaceTessellation);
  }

  buffer.updatePixels();
  buffer.endDraw();
  image(buffer, 0, 0); // draw our raster image on the screen
}

/*
 get array of triangles and draws them on the raster by
 calling draw2DTriangle()
 */
void drawTriangles(Triangle[] triangles) {
  for (Triangle t : triangles) {
    draw2DTriangle(t);  // Draw each triangle
  }
}

/*
Use the projected vertices to draw the 2D triangle on the raster.
 Several tasks need to be implemented:
 - cull degenerate or back-facing triangles
 - draw triangle edges using bresenhamLine()
 - draw normal vectors if needed
 - fill the interior using fillTriangle()
 */
void draw2DTriangle(Triangle t) {
  // Project the 3D vertices to 2D
  PVector[] projectedVertices = new PVector[3];
  for (int i = 0; i < 3; i++) {
    projectedVertices[i] = projectVertex(t.vertices[i]);
    if (projectedVertices[i] == null) {
      return; // Cull if any vertex projection fails
    }
  }

  // Cull back-facing triangles using the cross product in 2D
  if (!isFrontFacing(projectedVertices)) {
    return;
  }

  // Fill the triangle according to shading mode
  fillTriangle(t);  // function will handle the shading (flat, barycentric, etc.)

  // if outlines are enabled, draw the edges
  if (doOutline) {
    bresenhamLine((int) projectedVertices[0].x, (int) projectedVertices[0].y,
                  (int) projectedVertices[1].x, (int) projectedVertices[1].y);
    bresenhamLine((int) projectedVertices[1].x, (int) projectedVertices[1].y,
                  (int) projectedVertices[2].x, (int) projectedVertices[2].y);
    bresenhamLine((int) projectedVertices[2].x, (int) projectedVertices[2].y,
                  (int) projectedVertices[0].x, (int) projectedVertices[0].y);
  }

  // if normals are enabled, draw the normals
  if (doNormals) {
    drawNormals(t);
  }
}










/*
 Draw the normal vectors at each vertex and triangle center
 */
 
void drawNormals(Triangle t) {
  if (!doNormals) return;  // only draww if flag is set

  // Project the vertices to 2D space
  PVector[] projectedVertices = new PVector[3];
  for (int i = 0; i < 3; i++) {
    projectedVertices[i] = projectVertex(t.vertices[i]);
    if (projectedVertices[i] == null) {
      return; // Cull if any vertex projection fails
    }
  }

  // Calculate the center of the triangle
  PVector center = PVector.add(PVector.add(t.vertices[0], t.vertices[1]), t.vertices[2]).div(3);
  PVector projectedCenter = projectVertex(center);
  if (projectedCenter == null) {
    return;
  }

  // Average the normals
  PVector avgNormal = PVector.add(PVector.add(t.vertexNormals[0], t.vertexNormals[1]), t.vertexNormals[2]).div(3);
  avgNormal.mult(NORMAL_LENGTH);  // Scale normal for visibility
  PVector normalEnd = PVector.add(center, avgNormal);
  PVector projectedNormalEnd = projectVertex(normalEnd);

  if (projectedNormalEnd != null) {
    // make normal line at the center using bresenham
    bresenhamLine((int) projectedCenter.x, (int) projectedCenter.y,
                  (int) projectedNormalEnd.x, (int) projectedNormalEnd.y);
  }
}

/*
Fill the triangle according to the shading mode
*/
void fillTriangle(Triangle t) {
  if (shadingMode == ShadingMode.NONE) {
    return; // Do not fill the triangle
  }
  // Project the 3D vertices to 2D
  PVector[] projectedVertices = new PVector[3];
  for (int i = 0; i < 3; i++) {
    projectedVertices[i] = projectVertex(t.vertices[i]);
    if (projectedVertices[i] == null) {
      return; // Cull if any vertex projection fails
    }
  }

  // Compute the bounding box of the triangle
  int minX = (int) floor(min(projectedVertices[0].x, min(projectedVertices[1].x, projectedVertices[2].x)));
  int maxX = (int) ceil(max(projectedVertices[0].x, max(projectedVertices[1].x, projectedVertices[2].x)));
  int minY = (int) floor(min(projectedVertices[0].y, min(projectedVertices[1].y, projectedVertices[2].y)));
  int maxY = (int) ceil(max(projectedVertices[0].y, max(projectedVertices[1].y, projectedVertices[2].y)));

  // Loop over the bounding box
  for (int y = minY; y <= maxY; y++) {
    for (int x = minX; x <= maxX; x++) {
      // Perform point-in-triangle test
      if (pointInTriangle(x + 0.5f, y + 0.5f, projectedVertices)) {
        // Compute color according to shading mode
        float[] col = computeColorAtPixel(x + 0.5f, y + 0.5f, t, projectedVertices);
        setColor(col);
        setPixel(x + 0.5f, y + 0.5f);
      }
    }
  }
}

// Compute the color at a given pixel based on the current shading mode
float[] computeColorAtPixel(float x, float y, Triangle t, PVector[] projectedVertices) {
  float[] col = new float[3];
  switch (shadingMode) {
    case NONE:
      // No color
      break;
    case FLAT:
      col = FLAT_FILL_COLOR.clone();
      break;
    case BARYCENTRIC:
      float[] bary = computeBarycentricCoordinates(x, y, projectedVertices);
      col[0] = bary[0];
      col[1] = bary[1];
      col[2] = bary[2];
      break;
    case PHONG_FACE:
      // Compute Phong lighting at triangle center
      PVector center = PVector.add(PVector.add(t.vertices[0], t.vertices[1]), t.vertices[2]).div(3);
      PVector normal = computeFaceNormal(t);
      col = phong(center, normal, EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      break;
    case PHONG_VERTEX:
      // Compute Phong lighting at each vertex and average
      float[] c0 = phong(t.vertices[0], t.vertexNormals[0], EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      float[] c1 = phong(t.vertices[1], t.vertexNormals[1], EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      float[] c2 = phong(t.vertices[2], t.vertexNormals[2], EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      col[0] = (c0[0] + c1[0] + c2[0]) / 3;
      col[1] = (c0[1] + c1[1] + c2[1]) / 3;
      col[2] = (c0[2] + c1[2] + c2[2]) / 3;
      break;
    case PHONG_GOURAUD:
      // computate Phong lighting at each vertex, interpolate colors
      float[] baryGouraud = computeBarycentricCoordinates(x, y, projectedVertices);
      float[] c0G = phong(t.vertices[0], t.vertexNormals[0], EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      float[] c1G = phong(t.vertices[1], t.vertexNormals[1], EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      float[] c2G = phong(t.vertices[2], t.vertexNormals[2], EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      col[0] = baryGouraud[0]*c0G[0] + baryGouraud[1]*c1G[0] + baryGouraud[2]*c2G[0];
      col[1] = baryGouraud[0]*c0G[1] + baryGouraud[1]*c1G[1] + baryGouraud[2]*c2G[1];
      col[2] = baryGouraud[0]*c0G[2] + baryGouraud[1]*c1G[2] + baryGouraud[2]*c2G[2];
      break;
    case PHONG_SHADING:
      // calculate interpolated normal at pixel, then compute Phong lighting
      float[] baryPhongShading = computeBarycentricCoordinates(x, y, projectedVertices);
      PVector normalPS = PVector.mult(t.vertexNormals[0], baryPhongShading[0]);
      normalPS.add(PVector.mult(t.vertexNormals[1], baryPhongShading[1]));
      normalPS.add(PVector.mult(t.vertexNormals[2], baryPhongShading[2]));
      normalPS.normalize();

      // Interpolate position of shade
      PVector positionPS = PVector.mult(t.vertices[0], baryPhongShading[0]);
      positionPS.add(PVector.mult(t.vertices[1], baryPhongShading[1]));
      positionPS.add(PVector.mult(t.vertices[2], baryPhongShading[2]));

      col = phong(positionPS, normalPS, EYE, LIGHT, MATERIAL, PHONG_COLORS, PHONG_SHININESS);
      break;
  }
  // Clamp clamp clamp! color values to [0,1]
  col[0] = constrain(col[0], 0, 1);
  col[1] = constrain(col[1], 0, 1);
  col[2] = constrain(col[2], 0, 1);
  return col;
}

// Compute the face normal of the triangle
PVector computeFaceNormal(Triangle t) {
  PVector edge1 = PVector.sub(t.vertices[1], t.vertices[0]);
  PVector edge2 = PVector.sub(t.vertices[2], t.vertices[0]);
  PVector normal = edge1.cross(edge2);
  normal.normalize();
  return normal;
}

// Phong lighting model implementation
float[] phong(PVector p, PVector n, PVector eye, PVector light, float[] material, float[][] fillColor, float shininess) 
  {
  PVector L = PVector.sub(light, p).normalize(); // Light direction
  PVector V = PVector.sub(eye, p).normalize();   // View direction
  float NdotL = max(n.dot(L), 0);
  PVector R = PVector.sub(PVector.mult(n, 2 * NdotL), L).normalize(); // reflection direction
  float RdotV = max(R.dot(V), 0);
  float specularFactor = pow(RdotV, shininess);

  float[] result = new float[3];
  for (int i = 0; i < 3; i++) {
    float Ia = material[A] * fillColor[A][i];       // Ambient term
    float Id = material[D] * fillColor[D][i] * NdotL;  // Diffuse term
    float Is = material[S] * fillColor[S][i] * specularFactor; // Specular term
    result[i] = Ia + Id + Is;
  }
  return result;
}

// Sphere tessellation function i.e sphere to triangles
Triangle[] tessellateSphere(float r, int nPhi, int nTheta) {
    PVector[][] vertices = new PVector[nPhi + 1][nTheta];
    PVector[][] normals = new PVector[nPhi + 1][nTheta];

    float dPhi = PI / nPhi;
    float dTheta = TWO_PI / nTheta;

    // Generate vertices and normals
    for (int i = 0; i <= nPhi; i++) {
        float phi = i * dPhi;
        float sinPhi = sin(phi);
        float cosPhi = cos(phi);

        for (int j = 0; j < nTheta; j++) {
            float theta = j * dTheta;
            float sinTheta = sin(theta);
            float cosTheta = cos(theta);

            float x = r * sinPhi * cosTheta;
            float y = r * cosPhi;
            float z = r * sinPhi * sinTheta;

            vertices[i][j] = new PVector(x, y, z);

            // Normal vector
            float nx = sinPhi * cosTheta;
            float ny = cosPhi;
            float nz = sinPhi * sinTheta;

            normals[i][j] = new PVector(nx, ny, nz);
            normals[i][j].normalize();
        }
    }

    ArrayList<Triangle> triangles = new ArrayList<Triangle>();

    // Handle the north pole cap
    for (int j = 0; j < nTheta; j++) {
        PVector v0 = vertices[0][j]; // North pole
        PVector v1 = vertices[1][j];
        PVector v2 = vertices[1][(j + 1) % nTheta];

        PVector n0 = normals[0][j];
        PVector n1 = normals[1][j];
        PVector n2 = normals[1][(j + 1) % nTheta];

        // Ensure correct vertex order for CCW winding
        triangles.add(new Triangle(new PVector[] { v0, v2, v1 }, new PVector[] { n0, n2, n1 }));
    }

    // Handle the middle sections
    for (int i = 1; i < nPhi - 1; i++) {
        for (int j = 0; j < nTheta; j++) {
            PVector v0 = vertices[i][j];
            PVector v1 = vertices[i + 1][j];
            PVector v2 = vertices[i][(j + 1) % nTheta];
            PVector v3 = vertices[i + 1][(j + 1) % nTheta];

            PVector n0 = normals[i][j];
            PVector n1 = normals[i + 1][j];
            PVector n2 = normals[i][(j + 1) % nTheta];
            PVector n3 = normals[i + 1][(j + 1) % nTheta];

            // First triangle
            triangles.add(new Triangle(new PVector[] { v0, v2, v1 }, new PVector[] { n0, n2, n1 }));
            // Second triangle
            triangles.add(new Triangle(new PVector[] { v1, v2, v3 }, new PVector[] { n1, n2, n3 }));
        }
    }

    // Handle the south pole cap
    for (int j = 0; j < nTheta; j++) {
        PVector v0 = vertices[nPhi][j]; // South pole
        PVector v1 = vertices[nPhi - 1][j];
        PVector v2 = vertices[nPhi - 1][(j + 1) % nTheta];

        PVector n0 = normals[nPhi][j];
        PVector n1 = normals[nPhi - 1][j];
        PVector n2 = normals[nPhi - 1][(j + 1) % nTheta];

        // Ensure correct vertex order for CCW winding
        triangles.add(new Triangle(new PVector[] { v0, v1, v2 }, new PVector[] { n0, n1, n2 }));
    }

    return triangles.toArray(new Triangle[0]);
}
