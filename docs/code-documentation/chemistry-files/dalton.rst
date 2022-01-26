=========================================
Dalton files (``chemistry_files.dalton``)
=========================================

Files from the `Dalton software <http://daltonprogram.org/>`_.

.. note::

    About ``DALTON.PROP``:

    .. code-block:: text

        NORD = 0    energy (ground or excited)
               1    exp. value
               2    Linear response function
               3    Quadratic response function
               4    Cubic response function
              -1    ground - excited  transition matrix element, <0|x|i>
              -2    excited - excited transition matrix element, |<i|x|f>|
              -11   First order excited state property, <i|x|i>
              -20   <0|x|i><i|y|0>
              -21   w*<0|x|i><i|y|0>
              -22   (w_f - w_i)*|<i|x|f>|^2
              -30   D_pa
              -31   D_pe
              -32   D_pc
              -33   w1w2D_pa
              -34   w1w2D_pe
              -35   w1w2D_pc
               400  oscillator strength
               401  chemical shielding isotropic
               402  chemical shielding tensor

.. automodule:: qcip_tools.chemistry_files.dalton
    :members: