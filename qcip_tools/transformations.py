"""
Transformation (matrices) manipulation.

Constructed on top of `transforms3d <https://github.com/matthew-brett/transforms3d>`_ for the most part.
"""

import numpy


from transforms3d.euler import euler2mat
from transforms3d.axangles import axangle2mat


class TransformationMatrix:
    """4x4 matrix to represent a transformation
    """

    @staticmethod
    def identity():
        """Create identity 4x4 matrix

        :rtype: numpy.ndarray
        """
        return numpy.identity(4, dtype=numpy.float32)

    @staticmethod
    def translate(x, y, z):
        """Create a translation matrix

        :param x: translation x
        :type x: float
        :param y: translation y
        :param y: float
        :param z: translation z:
        :type z: float
        :rtype: numpy.ndarray
        """

        m = TransformationMatrix.identity()
        m[0:3, 3] = [x, y, z]
        return m

    @staticmethod
    def scale(x, y, z):
        """Create a translation matrix

        :param x: scale x
        :type x: float
        :param y: scale y
        :param y: float
        :param z: scale z
        :type z: float
        :rtype: numpy.ndarray
        """

        m = TransformationMatrix.identity()

        for i, j in zip([x, y, z], range(3)):
            m[j, j] = i

        return m

    @staticmethod
    def scale_uniform(factor):
        """Scale by a given factor

        :param factor: scale x
        :type factor: float
        :rtype: numpy.ndarray
        """

        return TransformationMatrix.scale(factor, factor, factor)

    @staticmethod
    def rotate_around_axis(axis, angle, is_normalized=False, in_degree=False):
        """Create a axis-angle rotation matrix

        :param axis: axis
        :type axis: numpy.ndarray
        :param angle: rotation angle (in degree)
        :type angle: float
        :param is_normalized: is the axis already normalized?
        :type is_normalized: bool
        :param in_degree: is the angle in degree ?
        :type in_degree: bool
        :rtype: numpy.ndarray
        """

        if in_degree:
            angle = numpy.radians(angle)

        m = TransformationMatrix.identity()
        m[:3, :3] = axangle2mat(axis, angle, is_normalized)
        return m

    @staticmethod
    def rotate_euler(psi, theta, chi, axes='sxyz', in_degree=False):
        """Rotate via Euler angle

        :param psi: first rotation angle
        :type psi: float
        :param theta: second rotation angle
        :type theta: float
        :param chi: third rotation angle
        :type chi: float
        :param axes: axes
        :type axes: str
        :param in_degree: is the angle in degree ?
        :type in_degree: bool
        :rtype: numpy.ndarray
        """

        if in_degree:
            psi = numpy.radians(psi)
            theta = numpy.radians(theta)
            chi = numpy.radians(chi)

        m = TransformationMatrix.identity()
        m[:3, :3] = euler2mat(psi, theta, chi, axes)
        return m


class ImmutableTransformable:
    """Transformation for non-mutable objects
    """

    def _apply_transformation(self, transformation):
        """Apply a transformation to an object

        :param transformation: the transformation
        :type transformation: numpy.ndarray
        :rtype: ImmutableTransformable
        """

        raise NotImplementedError()


class MutableTransformable:
    """Transformation for mutable objects
    """

    def _apply_transformation_self(self, transformation):
        """Apply a transformation to an object

        :param transformation: the transformation
        :type transformation: numpy.ndarray
        """

        raise NotImplementedError()


class ImmutableTranslatable(ImmutableTransformable):
    """Translatable objects"""

    def translate(self, x, y, z):
        """Translate object

        :param x: translation x
        :type x: float
        :param y: translation y
        :param y: float
        :param z: translation z:
        :type z: float
        :rtype: ImmutableTranslatable
        """

        return self._apply_transformation(TransformationMatrix.translate(x, y, z))


class MutableTranslatable(MutableTransformable):
    """Translatable objects"""

    def translate_self(self, x, y, z):
        """Translate object

        :param x: translation x
        :type x: float
        :param y: translation y
        :param y: float
        :param z: translation z:
        :type z: float
        """

        self._apply_transformation_self(TransformationMatrix.translate(x, y, z))


class ImmutableScalable(ImmutableTransformable):

    def scale(self, x, y, z):
        """Create a translation matrix

        :param x: scale x
        :type x: float
        :param y: scale y
        :param y: float
        :param z: scale z
        :type z: float
        :rtype: ImmutableScalable
        """

        return self._apply_transformation(TransformationMatrix.scale(x, y, z))

    def scale_uniform(self, factor):
        """Scale by a given factor

        :param factor: scale x
        :type factor: float
        :rtype: ImmutableScalable
        """

        return self._apply_transformation(TransformationMatrix.scale_uniform(factor))


class MutableScalable(MutableTransformable):

    def scale_self(self, x, y, z):
        """Create a translation matrix

        :param x: scale x
        :type x: float
        :param y: scale y
        :param y: float
        :param z: scale z
        :type z: float
        """

        self._apply_transformation_self(TransformationMatrix.scale(x, y, z))

    def scale_uniform_self(self, factor):
        """Scale by a given factor

        :param factor: scale x
        :type factor: float
        """

        self._apply_transformation_self(TransformationMatrix.scale_uniform(factor))


class ImmutableRotatable(ImmutableTransformable):

    def rotate_around_axis(self, axis, angle, is_normalized=False, in_degree=False):
        """Create a axis-angle rotation matrix

        :param axis: axis
        :type axis: numpy.ndarray
        :param angle: rotation angle (in degree)
        :type angle: float
        :param is_normalized: is the axis already normalized?
        :type is_normalized: bool
        :param in_degree: is the angle in degree ?
        :type in_degree: bool
        :rtype: ImmutableRotatable
        """

        return self._apply_transformation(
            TransformationMatrix.rotate_around_axis(axis, angle, is_normalized, in_degree))

    def rotate_euler(self, psi, theta, chi, axes='sxyz', in_degree=False):
        """Rotate via Euler angle

        :param psi: first rotation angle
        :type psi: float
        :param theta: second rotation angle
        :type theta: float
        :param chi: third rotation angle
        :type chi: float
        :param axes: axes
        :type axes: str
        :param in_degree: is the angle in degree ?
        :type in_degree: bool
        :rtype: ImmutableRotatable
        """

        return self._apply_transformation(TransformationMatrix.rotate_euler(psi, theta, chi, axes, in_degree))


class MutableRotatable(MutableTransformable):

    def rotate_around_axis_self(self, axis, angle, is_normalized=False, in_degree=False):
        """Create a axis-angle rotation matrix

        :param axis: axis
        :type axis: numpy.ndarray
        :param angle: rotation angle (in degree)
        :type angle: float
        :param is_normalized: is the axis already normalized?
        :type is_normalized: bool
        :param in_degree: is the angle in degree ?
        :type in_degree: bool
        """

        self._apply_transformation_self(
            TransformationMatrix.rotate_around_axis(axis, angle, is_normalized, in_degree))

    def rotate_euler_self(self, psi, theta, chi, axes='sxyz', in_degree=False):
        """Rotate via Euler angle

        :param psi: first rotation angle
        :type psi: float
        :param theta: second rotation angle
        :type theta: float
        :param chi: third rotation angle
        :type chi: float
        :param axes: axes
        :type axes: str
        :param in_degree: is the angle in degree ?
        :type in_degree: bool
        """

        self._apply_transformation_self(TransformationMatrix.rotate_euler(psi, theta, chi, axes, in_degree))
