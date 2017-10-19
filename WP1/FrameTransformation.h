#ifndef FRAMETRANSFORMATION_H
#define FRAMETRANSFORMATION_H

Eigen::Vector6d FrameTransformationSE(long double x, long double y, long double xdot, long double ydot);

Eigen::Vector6d FrameTransformationEM(long double x, long double y, long double xdot, long double ydot);

#endif // FRAMETRANSFORMATION_H
