#ifndef ADD_H
#define ADD_H

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

//OpenGL stuff
GLfloat angle = 0.0;
void sphere(void) {
	glColor3f(1.0,0.0,0.0);
	glRotatef(angle,1.0,0.0,0.0);
	glutSolidSphere(1.0,20,20);
}


void display(void) {
	glClearColor(0.0,0.0,0.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	//gluLookAt(0.0,0.0,5.0,0.0,0.0,0.0,0.0,1.0,0.0);
	sphere();
	//angle++;
	glFlush();
}
#endif /* ADD_H */
