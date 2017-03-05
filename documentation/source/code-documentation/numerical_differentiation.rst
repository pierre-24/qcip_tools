====================================================================
Numerical differentiation (``qcip_tools.numerical_differentiation``)
====================================================================


Mathematical background
-----------------------

Numerical derivation
++++++++++++++++++++

Theory
______

This an application of the `Neville's algorithm <https://en.wikipedia.org/wiki/Neville%27s_algorithm>`_.


Let :math:`M(x)` be a `Maclaurin series expansion <https://en.wikipedia.org/wiki/Taylor_series>`_ of a function :math:`f(x)`,

.. math::
    :label: expansion

    M(x) = \sum_ {n=0}^{n_{max}} \frac{A_n}{n!} x^n,

where :math:`A_n=f^{(n)}(0)`, the value of the :math:`n`-order derivative of :math:`f(x)` at :math:`x=0`. If the expansion is not truncated, :math:`n_{max}=\infty`.

Given a small value :math:`h>0` and :math:`I_{a,q} \in \mathbb{R}`,

.. math::

    M(h\,I_{a,q}) = \sum_{n=0}^{n_{max}} I_{a,q}^n\,\frac{A_n}{n!}\,h^n

where

.. math::

    \forall q \in \mathbb{N}: I_{a,q} = \left\{\begin{matrix}-a^{|q|-1} &\text{if }q< 0,\\0 &\text{if } q= 0, \\ a^{q-1} & \text{if }q> 0.\end{matrix}\right.

Therefore, :math:`\forall h \in \mathbb{R}^+:\forall k\in\mathbb{N}^0:h\,I_{a, k+1}=a^k\,h` defines a geometric progression where :math:`h` is the smallest value and :math:`a` is the common ratio (or scale factor).

An estimate of an :math:`A_d` factor [value of :math:`f^{(d)}(0)`] determined by using :math:`h`, is given by

.. math::

    A_d(h) = \frac{d!}{h^d}\,M^{(d)}(h) + \mathcal{O}(x^{d+1}),

where :math:`M^{(d)}(h)` is the value at :math:`h` of the :math:`d`-order derivative of the Maclaurin series.

Let's assume that one can determine :math:`A_d(h)` with a given error :math:`\mathcal{O}(h^{d+p}), p\in\mathbb{N}^+` by computing

.. math::

    \begin{align}
    A_d(h) &= \frac{d!}{h^d}\,\sum_{q=q_{min}}^{q_{max}} C_q\,M(h\,I_{a,q})+  \mathcal{O}(h^{d+p}), \\
        &= \frac{d!}{h^d}\,\sum_{q=q_{min}}^{q_{max}} C_q\, \sum_{n=0}^{n_{max}} I_{a,q}^n\,\frac{A_n}{n!}\,h^n+  \mathcal{O}(h^{d+p}), \\
        &= \frac{d!}{h^d}\,\sum_{n=0}^{n_ {max}} \left[\sum_{q=q_{min}}^{q_{max}} I_{a,q}^n\,C_q\right]\,\frac{A_n}{n!}\,h^n +  \mathcal{O}(h^{d+p}).
    \end{align}

In order for the equality to be satisfied, it is necessary that:

.. math::
    :label: sys

    \sum_{q=q_{min}}^{q_{max}} I_{a,q}^n\,C_q = \left\{
        \begin{array}{ll}
            1 & \text{if }n=d,\\
            0 & \text{otherwise.}
        \end{array}\right.

This defines a set of :math:`n_{max}` linear equations with :math:`q_{max}-q_{min}+1` unknowns (:math:`C_q`).
In order to get a unique solution, :math:`q_{max}-q_{min}+1\geqslant d` and therefore, one can set :math:`n_{max}=q_{max}-q_{min}+1= d+p`.
:math:`q_{min}` and :math:`q_{max}` take different values, depending of the approximation:

+------------------+----------------------------------------------------+--------------------------------------------------+
| Type             | :math:`i_{min}`                                    | :math:`i_{max}`                                  |
+==================+====================================================+==================================================+
| Forward (``F``)  | 0                                                  | :math:`d+p-1`                                    |
+------------------+----------------------------------------------------+--------------------------------------------------+
| Backward (``B``) | :math:`-d+p-1`                                     | 0                                                |
+------------------+----------------------------------------------------+--------------------------------------------------+
| Centered (``C``) | :math:`-\left\lfloor\frac{d+p-1}{2}\right\rfloor`  | :math:`\left\lfloor\frac{d+p-1}{2}\right\rfloor` |
+------------------+----------------------------------------------------+--------------------------------------------------+

Note that :math:`p` is always even in the case of a centered derivative approximation.

The system of linear equation given by Eq. :eq:`sys` is solved. Once the :math:`C_q`'s are determined, :math:`A_d(h)` is obtained with a precision :math:`\mathcal{O}(h^{d+p})` by computing

.. math::
    :label: ad

    A_d(h) = \frac{d!}{h^d}\,\sum_{q=q_{min}}^{q_{max}} C_q\,M(h\,I_{a,q}).

It can be simply showed that the coefficients for a multivariate function is a tensor product of the coefficients for the univariates approximations. In this case :math:`A_d(\mathbf{h})` is a symmetric tensor, while :math:`\mathbf{h}` is a vector.

.. note::

    In pratice, it wouldn't change anything if a Taylor series (not centered at :math:`x=0`) is used instead. It also holds for power series if one discard the :math:`d!` in :eq:`ad`.

.. warning::

    If values of :math:`M(h)` are known with a certain precision :math:`\delta y`, the absolute error on :math:`A_d(h)` cannot be lower than :math:`\delta y \times h^{-d}`, no matter the order of magnitude of :math:`A_d(h)`.
    If there is no physical (or chemical) limitation to the precision, :math:`\delta y` is then determined by the precision of the representation,
    see the `python documentation on floating point arithmetic <https://docs.python.org/3.6/tutorial/floatingpoint.html>`_.


Implementation
______________

The :math:`I_{a,q}` function is named ``ak_shifted()`` in the python code.

In the code, to works with the romberg triangles defined below, :math:`h=a^{k}h_0`, where :math:`k` is the minimal amplitude, while :math:`h_0` is the minimal value.

The ``Coefficients`` class takes care of the computation of the different coefficients, storing them in ``Coefficients.mat_coefs``.
The :math:`\frac{d!}{h^d}` part can be computed via the ``Coefficients.prefactor()`` function.

The ``compute_derivative_of_function()`` function does the whole computation (for uni and multivariate functions) and gives :math:`A_d(h)` (which is the component of a tensor in the case of multivariate functions).
A ``scalar_function(fields, h_0, **kwargs)`` function must be defined, which must return the (scalar) value of the (multivariate) Maclaurin expansion :math:`M()` given ``fields``, which are the different :math:`q` for each variable.

.. warning::

    In the code, :math:`k` is added (or substracted) to :math:`q` so that ``scalar_function()`` receives the :math:`q` values with respect to :math:`h_0` and is expected to returns the value of :math:`M(I_{a,q}h_0)`.



Romberg's scheme (Richardson extrapolation)
+++++++++++++++++++++++++++++++++++++++++++

Theory
______


In quantum chemistry, the `Richardson extrapolation <https://en.wikipedia.org/wiki/Richardson_extrapolation>`_ is also known as the Romberg differentiation scheme.
Like its integration counterpart, it allows to refine the value of a numerical derivatives by computing multiple estimate of the derivates.

Given Eq. :eq:`expansion` and a small :math:`h \in \mathbb{R}^+`, :math:`r, m \in \mathbb{N}^+`, :math:`k \in \mathbb{N}^0`,

.. math ::

    \begin{align}
        a^{rm}\,M(a^k h) - M(a^{k+1} h) &= a^{rm}\,\sum_ {n=0}^{n_{max}} a^{nk}\,\frac{A_n}{n!} h^n - \sum_ {n=0}^{n_{max}} a^{nk+n}\,\frac{A_n}{n!} h^n \\
            &= \sum_ {n=0}^{n_{max}} [a^{nk+rm}-a^{nk+n}]\frac{A_n}{n!} h^n.
    \end{align}



It is therefore possible to remove the :math:`s`-power term from the expansion by choosing :math:`s=rm`.

Let :math:`H_{k,0}\equiv A_d(h=a^k h_0)` with :math:`h_0` the minimal value, the following recurrence relation can be defined:

.. math::
    :label: romberg

    H_{k,m+1} = \frac{a^{rm}\,H_{k,m}-H_{k+1,m}}{a^{rm}-1} + \mathcal{O}(h^{rm}),


where :math:`m` is the number of *refinement* steps or iterations. To achieve such a :math:`\mathcal{O}(h^{rm})` precision, it is required to know :math:`m+1` values of :math:`H_{k,0}` with :math:`k\in [0;m+1]`.
In general, :math:`r` should be equal to 1 to remove the :math:`m`-power contamination at iteration :math:`m`, but in the case of centered derivatives, every odd-power term vanishes, so one can use :math:`r=2` to achieve a faster convergence by removing the :math:`2m`-power contamination at iteration :math:`m`.

A romberg triangle is obtained, with the following shape:

.. code-block:: text

    k=  | h=    |  m=0     m=1     m=2    ...
    ----+-------+-----------------------------
    0   | a⁰*h0 |  val01   val11   val20  ...
    1   | a¹*h0 |  val02   val12   ...
    2   | a²*h0 |  val03   ...
    ... | ...   |  ...

If the minimal value :math:`h_0` is chosen well, the value of the derivative should be the rightmost value. In practice, there is a value window with lower bound (for which there are too large round-off errors) and larger bound (when higher order terms makes the number of :math:`a^kh_0` required large, along with quantum chemistry related reasons which lower that value).
Without knowledge of this ideal :math:`h_0` window, it is necessary to carry out an analysis of the triangle to select the "best" value. Two quantities are usefull:

.. math::

    \begin{align}
        &\varepsilon_k(m) = H_{k+1,m} - H_{k,m}, \\
        &\varepsilon_{m}(k) = H_{k,m+1} - H_{k,m},
    \end{align}

where :math:`\varepsilon_k(m)` is the amplitude error at a given iteration :math:`m` and :math:`\varepsilon_m(k)` is the iteration error for a given amplitude :math:`k`.

This "best" value is chosen according to the following flowchart:

.. figure:: /_static/flowchart_romberg_selection.png
    :align: center

    Flowchart to select the "best" value in a Romberg triangle, adapted from the text in M. de Wergifosse *et al*. *Int. J. Quant. Chem.* **114**, 900 (2014).

Note that the original paper does not give advices on what the threshold value should be, but it clearly depends on the order of magnitude of :math:`A_d(h)` and the available (and/or required) precision.

.. note::

    In the case where :math:`A_d(\mathbf{h}=a^k\mathbf{h}_0)` is a tensor (because the corresponding function is multivariate), a Romberg triangle must be carried out for each component of this tensor.
    But since the vector is symmetric, it can reduce the number of components to carry out.

Implementation
______________

The ``RombergTriangle`` class implement the Romberg's triangle algorithm. The constructors expect the value corresponding to increasing ``k`` values to be passed as argument.

The :math:`H_{k,m}` are available as ``RombergTriangle.romberg_triangle[k, m]``. Note that the values for which :math:`k \geqslant k_{max} - m` are set to zero.

The ``find_best_value()`` implement the flowchart showed above. To help the (eventual) debugging and understanding, a ``verbose`` option is available.

For example, for an univariate function :math:`F(x)`, compute :math:`\frac{dF(x)}{dx}`:

.. code-block:: python

    a = 2.0
    # Forward first order derivative:
    c = numerical_differentiation.Coefficients(1, 1, ratio=a, method='F')
    derivative_values = []

    # compute dF(x) / dx
    for k in range(5):
        derivative_values.append(numerical_differentiation.compute_derivative_of_function(
            [(c, 0)],  # derivative with respect to the first (and only) variable
            univariate_function,  # assume defined
            k,
            0.001,  # f0
            1  # differentiation space size is 1 because univariate function
        ))

    t = numerical_differentiation.RombergTriangle(derivative_values, ratio=a)
    position, value, iteration_error = t.find_best_value(threshold=1e-5, verbose=True)

For a bivariate function :math:`F(x, y)`, compute :math:`\frac{\partial^3 F(x, y)}{\partial x \partial y^2}`:

.. code-block:: python

    a = 2.0
    # centered first and second order derivative:
    c1 = numerical_differentiation.Coefficients(1, 2, ratio=a, method='C')
    c2 = numerical_differentiation.Coefficients(2, 2, ratio=a, method='C')
    derivative_values = []

    # compute d³F(x,y) / dxdy²
    for k in range(5):
        derivative_values.append(numerical_differentiation.compute_derivative_of_function(
            [(c1, 0), (c2, 1)],  # first and second order derivative with respect to the first and second variable
            bivariate_function,  # assume defined
            k,
            0.001,  # f0
            2  # differentiation space size is 2
        ))

    t = numerical_differentiation.RombergTriangle(derivative_values, ratio=a, r=2)  # r=2 because centered
    position, value, iteration_error = t.find_best_value()



Sources
-------

+ D\. Eberly, « Derivative Approximation by Finite Difference ». From the `Geometric Tools documentation <https://www.geometrictools.com/Documentation/Documentation.html>`_ (last consultation: March 4, 2017).
+ J.N. Lyness *et al*. *Numer. Math.* **8**, 458 (1966).
+ L.F. Richardson *et al*. *Philosophical Transactions of the Royal Society of London. Series A, Containing Papers of a Mathematical or Physical Character* **226**, 299 (1927).
+ M\. de Wergifosse *et al*. *Int. J. Quant. Chem.* **114**, 900 (2014).
+ A.A.K. Mohammed *et al*. *J. Comp. Chem.* **34**, 1497 (2013).

API documentation
-----------------

.. automodule:: qcip_tools.numerical_differentiation
    :members: