/*
Use Bresenham's line algorithm to draw a line on the raster between
 the two given points. Modify the raster using setPixel() ONLY.
 */
void bresenhamLine(int fromX, int fromY, int toX, int toY) {
    
    int dx = abs(toX - fromX);
    int dy = abs(toY - fromY);
    
    int sx = fromX < toX ? 1 : -1;  // Step direction for x
    int sy = fromY < toY ? 1 : -1;  // Step direction for y

    int err = dx - dy;  // Initial error

    while (true) {
        setPixel(fromX, fromY);  // Set the current pixel to the desired color
        
        if (fromX == toX && fromY == toY) // Exit loop when the end point is reached
          break;  
        
        int e2 = 2 * err;
        
        if (e2 > -dy) {
            err -= dy;
            fromX += sx;  // Move in x direction
        }
        if (e2 < dx) {
            err += dx;
            fromY += sy;  // Move in y direction
        }
    }
}

/*
Don't change anything below here
 */

final int LENGTH_X = 125;
final int LENGTH_Y = 125;
final int LENGTH_DIAGONAL = 52;

void testBresenham() {
  final color WHITE = color(1f, 1f, 1f);
  final color RED = color(1f, 0f, 0f);

  final int CENTER_OFFSET_X = 125;
  final int CENTER_OFFSET_Y = 125;

  buffer.updatePixels(); // display everything drawn so far

  buffer.stroke(RED);
  ComparisonLines comp = new ComparisonLines();
  comp.drawAllQuadrants(CENTER_OFFSET_X, CENTER_OFFSET_Y);

  buffer.loadPixels(); // switch back to editing the raster
  setColor(WHITE);

  // use the implementation of Bresenham's algorithm
  BresenhamLines bres = new BresenhamLines();
  bres.drawAllQuadrants(CENTER_OFFSET_X, CENTER_OFFSET_Y);
}

abstract class TestPattern {
  void drawAllQuadrants(int centerOffsetX, int centerOffsetY) {
    for (int signX=-1; signX<=1; signX+=2) {
      int centerX = signX*centerOffsetX;
      for (int signY=-1; signY<=1; signY+=2) {
        int centerY = signY*centerOffsetY;
        drawPattern(centerX, centerY);
      }
    }
  }

  void drawPattern(int centerX, int centerY) {
    drawAxes(centerX, centerY);
    drawDiagonals(centerX, centerY);
  }

  void drawAxes(int startX, int startY) {
    for (int sign=-1; sign<=1; sign+=2) {
      int endXHorizontal = startX + sign*LENGTH_X;
      int endYVertical = startY + sign*LENGTH_Y;
      drawLine(startX, startY, endXHorizontal, startY);
      drawLine(startX, startY, startX, endYVertical);
    }
  }

  void drawDiagonals(int startX, int startY) {
    for (int signX=-1; signX<=1; signX+=2) {
      int endXHorizontal = startX + signX*LENGTH_X;
      int endXDiagonal = startX + signX*LENGTH_DIAGONAL;
      for (int signY=-1; signY<=1; signY+=2) {
        int endYVertical = startY + signY*LENGTH_Y;
        int endYDiagonal = startY + signY*LENGTH_DIAGONAL;
        drawLine(startX, startY, endXDiagonal, endYVertical);
        drawLine(startX, startY, endXHorizontal, endYDiagonal);
      }
    }
  }

  abstract void drawLine(int startX, int startY, int endX, int endY);
}

class ComparisonLines extends TestPattern {
  void drawLine(int startX, int startY, int endX, int endY) {
    final int SMALL_SHIFT = 3;

    // shift left/right or up/down
    int xDir = -Integer.signum(endY - startY);
    int yDir = Integer.signum(endX - startX);

    int pStartX = rasterToProcessingX(startX + xDir*SMALL_SHIFT);
    int pStartY = rasterToProcessingY(startY + yDir*SMALL_SHIFT);

    int pEndX = rasterToProcessingX(endX + xDir*SMALL_SHIFT);
    int pEndY = rasterToProcessingY(endY + yDir*SMALL_SHIFT);

    buffer.line(pStartX, pStartY, pEndX, pEndY);
  }

  int rasterToProcessingX(int rx) {
    return width/2 + rx;
  }

  int rasterToProcessingY(int ry) {
    return height/2 - ry;
  }
}

class BresenhamLines extends TestPattern {
  void drawLine(int startX, int startY, int endX, int endY) {
    bresenhamLine(startX, startY, endX, endY);
  }
}
