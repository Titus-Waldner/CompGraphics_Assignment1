Assignment 1 - Computer Graphics
Overview

This project implements fundamental computer graphics algorithms using Processing 4.3. The assignment focuses on drawing and rendering 3D objects pixel by pixel, without using built-in drawing functions like point(), line(), or triangle().
Features Implemented
Bresenham’s Line Algorithm

    Implemented to draw lines in all eight octants using the setPixel() function.
    Validated by comparing results against Processing’s line() function.
    
![image](https://github.com/user-attachments/assets/d0a26bd3-4a38-4f57-86f1-c31a10b17a10)

3D Triangle Rendering

    Created and transformed a single 3D triangle in space.
    Used perspective projection with the provided projectVertex() function.
    Implemented back-face culling to remove degenerate or hidden triangles.
    Visualized normal vectors for debugging.
    

Surface Tessellation

    Constructed a tessellated 3D shape using triangles.
    Supported multiple shapes including a hyperboloid, sphere, and cone.
    Allowed adjustable tessellation parameters for different levels of detail.

![image](https://github.com/user-attachments/assets/2bbaa611-1b26-4b94-8291-531fd8cec2f3)


Scan-line Triangle Fill Algorithm

    Implemented scan-line filling to render solid triangles.
    Supported flat shading and barycentric coordinate visualization.

Lighting and Shading

    Implemented Phong lighting calculations for rendering.
    Supported face-level Phong lighting, vertex-level Phong lighting, Gouraud shading, and Phong shading.

![image](https://github.com/user-attachments/assets/3f17ccf1-6c90-4e62-819a-ad2feeacbe88)


Code Structure and Optimization

    Used structured classes for geometric objects.
    Applied mathematical operations with Processing’s PVector and custom utility functions.
    Followed modular programming practices to keep functions organized and readable.

Usage

To run the project, open the .pde files in Processing 4.3 and execute the main script. The program responds to predefined hotkeys for different drawing and shading modes.
Submission

All necessary .pde files are included in the repository. The implementation follows the provided template structure without modifying hotkeys or core functionality.
