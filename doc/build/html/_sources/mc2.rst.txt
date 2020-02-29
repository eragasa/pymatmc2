Multi-Cell Monte Carlo (MC)^2
=============================

I am currently working on pymatmc2 which implements the methodology of Niu, Rao, Windl, and Ghazisaeidi :cite:`niu2019_mc2_flip`

.. math::
   :label: lever_rule
   
   x_1^1 f_1 + x_2^1 f_2 = c^1
   
   x_2^2 f_2 + x_2^2 f_2 = c^2

In linear algebra form, Equation :eq:`lever_rule` becomes :math:`\bm{X}\bm{f}=\bm{c}`, so that :math:`\bm{f}=\bm{X}^{-1}\bm{f}` determines the molar ratios.  This approach preserves the initial stoichiometry.   

The probability of acceptance is the same as an isothemal-isobaric ensemble with the acceptane criteria for a move from :math:`k` to :math:`k+1`

.. math::
   :label: p_accept

   p_{\mathrm{accept}} = min\left[
       1, \exp\left(\beta \Delta H + N \sum_{i=1}^m \ln \frac{V_i(k+1)}{V_i(k)}\right)
   \right]
.. bibliography:: references.bib