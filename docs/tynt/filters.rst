
Available Filters
-----------------

The table below lists the filters available by default by calling :py:meth:`~tynt.FilterGenerator.reconstruct`
using the included, light-weight filter archive in ``fft.fits``. All filters provided by the `SVO Filter Profile Service
<http://svo2.cab.inta-csic.es/theory/fps/>`_ are accessible via ``tynt`` with the
:py:meth:`~tynt.FilterGenerator.download_true_transmittance` method, which will access the internet and download
the transmittance curve at its "native" resolution, without the FFT approximation.

If the filters you need are not listed below but you'd like offline access to the ``tynt`` interface, you can
always create your own archive locally, with a custom set of filters, including all available filters. For a
tutorial on downloading your own custom archive, see :doc:`customarchive`.


To access a filter from the table below, combine the columns like this: ``col1/col2.col3``. For example, if you want
the 2MASS J filter in the first row below, you would call the following in ``tynt``:

.. code-block:: python

    from tynt import FilterGenerator

    f = FilterGenerator()
    filt = f.reconstruct('2MASS/2MASS.J')

.. raw:: html

    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th></th>
          <th>Filters</th>
        </tr>
        <tr>
          <th></th>
          <th>Set</th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>2MASS</th>
          <th>2MASS</th>
          <td>J, H, Ks</td>
        </tr>
        <tr>
          <th>CAHA</th>
          <th>BUSCA</th>
          <td>u, v, b, y</td>
        </tr>
        <tr>
          <th>CHEOPS</th>
          <th>CHEOPS</th>
          <td>band</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">CTIO</th>
          <th>SOI</th>
          <td>bessel_U, bessel_B, bessel_V, bessel_R, bessel_I, strom_u, strom_v, strom_b, strom_y</td>
        </tr>
        <tr>
          <th>Y4KCam</th>
          <td>Rc, Ic</td>
        </tr>
        <tr>
          <th rowspan="5" valign="top">GAIA</th>
          <th>Gaia2m</th>
          <td>Gbp_faint, Gbp_bright, G, Grp</td>
        </tr>
        <tr>
          <th>GAIA2r</th>
          <td>Gbp, G, Grp</td>
        </tr>
        <tr>
          <th>GAIA2</th>
          <td>Gbp, G, Grp</td>
        </tr>
        <tr>
          <th>GAIA3</th>
          <td>Gbp, G, Grp, Grvs</td>
        </tr>
        <tr>
          <th>GAIA0</th>
          <td>Gbp, G, Grp</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">GCPD</th>
          <th>Johnson</th>
          <td>U_Landolt, U_ADPS, U_Mendoza, U, U_Straizys, B_Landolt, B, B_ADPS, B_Mendoza, B_Straizys, V_Landolt, V_Straizys, V, V_Mendoza, V_ADPS, R_Landolt, R_Mendoza, R, I_Landolt, I_Mendoza</td>
        </tr>
        <tr>
          <th>JHKLMN</th>
          <td>U, B, V, R, I</td>
        </tr>
        <tr>
          <th>Cape</th>
          <td>Uc, Uc_ADPS, B_ADPS, V_ADPS</td>
        </tr>
        <tr>
          <th>Eggen</th>
          <td>R, I</td>
        </tr>
        <tr>
          <th>Cousins</th>
          <td>R_ADPS, R, I_ADPS, I</td>
        </tr>
        <tr>
          <th>Stromgren</th>
          <td>u_ADPS, u, v_ADPS, v, b_ADPS, b, y, y_ADPS</td>
        </tr>
        <tr>
          <th rowspan="3" valign="top">Generic</th>
          <th>Johnson</th>
          <td>U, B, V, R, I, J, M</td>
        </tr>
        <tr>
          <th>Cousins</th>
          <td>R, I</td>
        </tr>
        <tr>
          <th>Stromgren</th>
          <td>u, v, b, y</td>
        </tr>
        <tr>
          <th rowspan="25" valign="top">HST</th>
          <th>STIS_FUV</th>
          <td>F25LYA, F25LYA_G140L, F25LYA_G140M, 25MAMA_G140M, F25ND3_G140M, 25MAMA_G140L, F25ND3_G140L, F25ND5_G140M, F25ND5_G140L, F25NDQ2_G140M, 25MAMA, F25ND3, F25NDQ3_G140M, F25NDQ2_G140L, F25NDQ3_G140L, F25NDQ1_G140M, F25ND5, F25NDQ1_G140L, F25NDQ4_G140M, F25NDQ4_G140L, F25NDQ2, F25NDQ3, F25NDQ1, F25SRF2_G140M, F25SRF2_G140L, F25NDQ4, F25SRF2, F25QTZ_G140L, F25QTZ_G140M, F25QTZ</td>
        </tr>
        <tr>
          <th>ACS_SBC</th>
          <td>F122M, F115LP, PR110L, F125LP, PR130L, F140LP, F150LP, F165LP</td>
        </tr>
        <tr>
          <th>WFPC1-WF</th>
          <td>G200M2, F194W, F230W, F284W, F336W, F368M, F375N, F157W, F413M, F437N, F439W, F469N, F487N, F492M, F502N, G450, F517N, F122M, F547M, G200, F555W, F569W, F588N, F606W, F622W, F631N, F648M, F656N, F658N, F664N, F673N, F675W, F8ND, POL0, POL120, POL60, F128LP, F702W, F718M, G800, F791W, F814W, F875M, F889N, F725LP, F785LP, F850LP, F1083N, F1042M</td>
        </tr>
        <tr>
          <th>WFPC1-PC</th>
          <td>G200M2, F194W, F230W, F284W, F336W, F368M, F375N, F413M, F437N, F439W, F469N, F487N, F492M, F502N, G450, F517N, F157W, F547M, F555W, F569W, F588N, G200, F606W, F622W, F122M, F631N, F648M, F656N, F658N, F664N, F673N, F675W, F702W, F718M, F8ND, POL0, POL120, POL60, F128LP, G800, F791W, F814W, F875M, F889N, F725LP, F785LP, F850LP, F1083N, F1042M</td>
        </tr>
        <tr>
          <th>WFPC2-WF</th>
          <td>F160BW, F185W, F170W, F218W, F255W, F157W, F300W, F336W, F343N, F122M, F375N, FQUVN33, FQUVN_B, F390N, FQUVN_C, FQUVN_D, F380W, F410M, F439W, F437N, F450W, F467M, F469N, F487N, F502N, FQCH4N_D, F547M, F555W, F569W, F588N, F606W, FQCH4N33, F622W, F631N, F656N, F658N, F130LP, F165LP, F673N, F675W, F702W, POLQ, POLQ_90, POLQ_45, FQCH4N_B, F791W, F953N, F814W, F785LP, FQCH4N_C, F850LP, F1042M</td>
        </tr>
        <tr>
          <th>WFPC2-PC</th>
          <td>F160BW, F185W, F170W, F218W, F255W, F157W, F300W, F336W, F343N, F122M, F375N, FQUVN, F390N, F380W, F410M, F439W, F437N, F450W, F467M, F469N, F487N, F502N, FQCH4N, F547M, F555W, F569W, F588N, F606W, F622W, F631N, F656N, F658N, F130LP, F165LP, F673N, F675W, F702W, POLQ, F791W, F953N, F814W, F785LP, F850LP, F1042M</td>
        </tr>
        <tr>
          <th>HSP_UV1</th>
          <td>F145M_A, F145M_B, F135W_A, F152M_A, F152M_B, F135W_B, PRISM_BLUE, F122M_A, F122M_B, F184W_A, F184W_B, F218M_A, F218M_B, F220W_A, F220W_B, F240W_A, F240W_B, F140LP_A, F140LP_B, F248M_A, F248M_B, PRISM_RED, F278N_A, F278N_B</td>
        </tr>
        <tr>
          <th>HSP_UV2</th>
          <td>F145M_A, F145M_B, PRISM_BLUE, F152M_A, F152M_B, F122M_A, F122M_B, F184W_A, F184W_B, F179M_A, F179M_B, F218M_A, F218M_B, F140LP_A, F140LP_B, F160LP_A, F160LP_B, F248M_A, F248M_B, F262M_A, F262M_B, PRISM_RED, F284M_A, F284M_B, F278N_A, F278N_B</td>
        </tr>
        <tr>
          <th>FOC_F48</th>
          <td>F140W, F150W, F175W, F220W, F195W, F275W, F342W, F430W, PRISM1, PRISM3, PRISM2, F130LP, F180LP, F305LP</td>
        </tr>
        <tr>
          <th>FOC_F96</th>
          <td>F140W, F130M, F170M, F175W, F210M, F190M, F120M, F152M, F220W, F231M, F165W, F140M, F253M, F278M, F275W, F307M, F320W, F342W, F195W, F346M, F372M, F410M, F430W, F437M, F470M, F486N, F502M, F501N, F550M, F600M, F6ND, F2ND, F1ND, PRISM1, PRISM2, F130LP, POL0, POL120, POL60, F4ND, F8ND, F370LP, F480LP, F630M</td>
        </tr>
        <tr>
          <th>HSP_VIS</th>
          <td>F184W_A, F184W_B, PRISM_BLUE, F240W_A, F240W_B, F262M_A, F262M_B, F355M_A, F355M_B, F419N_A, F419N_B, F450W_A, F450W_B, F160LP_A, F160LP_B, F551W_A, F551W_B, PRISM_RED, F400LP_A, F400LP_B, F620W_A, F620W_B</td>
        </tr>
        <tr>
          <th>STIS_NUV</th>
          <td>F25CIII_G230L, F25CIII_G230M, F25CIII, F25CIII_PRISM, F25CN182, F25CN182_PRISM, F25CN182_G230L, F25CN182_G230M, 25MAMA, F25SRF2, 25MAMA_PRISM, F25NDQ1, F25QTZ, F25SRF2_PRISM, F25NDQ1_PRISM, 25MAMA_G230L, F25NDQ2, F25QTZ_PRISM, F25QTZ_G230L, F25SRF2_G230L, 25MAMA_G230M, F25QTZ_G230M, F25SRF2_G230M, F25ND3, F25NDQ1_G230L, F25NDQ1_G230M, F25NDQ2_PRISM, F25ND3_PRISM, F25NDQ2_G230L, F25NDQ2_G230M, F25NDQ3, F25ND3_G230L, F25ND3_G230M, F25NDQ3_PRISM, F25NDQ3_G230L, F25NDQ3_G230M, F25NDQ4, F25NDQ4_G230L, F25NDQ4_G230M, F25NDQ4_PRISM, F25CN270_G230L, F25CN270, F25CN270_PRISM, F25CN270_G230M, F25ND5_G230M, F25ND5_G230L, F25ND5, F25MGII, F25MGII_PRISM, F25MGII_G230L, F25MGII_G230M, F25ND5_PRISM</td>
        </tr>
        <tr>
          <th>WFC3_UVIS2</th>
          <td>F218W, FQ232N, F225W, FQ243N, F275W, F280N, F300X, F336W, F343N, F373N, FQ378N, FQ387N, F390M, F390W, F395N, F410M, FQ422M, F438W, FQ436N, FQ437N, G280, F467M, F469N, F475W, F487N, FQ492N, F502N, F475X, FQ508N, F555W, F547M, FQ575N, F606W, F200LP, FQ619N, F621M, F625W, F631N, FQ634N, F645N, F350LP, F656N, F657N, F658N, F665N, FQ672N, FQ674N, F673N, F680N, F689M, FQ727N, FQ750N, F763M, F600LP, F775W, F814W, F845M, FQ889N, FQ906N, F850LP, FQ924N, FQ937N, F953N</td>
        </tr>
        <tr>
          <th>WFC3_UVIS1</th>
          <td>F218W, FQ232N, F225W, FQ243N, F275W, F280N, F300X, F336W, F343N, F373N, FQ378N, FQ387N, F390M, F395N, F390W, F410M, FQ422M, F438W, FQ436N, FQ437N, F467M, F469N, G280, F475W, F487N, FQ492N, F502N, F475X, FQ508N, F555W, F547M, FQ575N, F606W, F200LP, FQ619N, F621M, F625W, F631N, FQ634N, F645N, F350LP, F656N, F657N, F658N, F665N, FQ672N, FQ674N, F673N, F680N, F689M, FQ727N, FQ750N, F763M, F600LP, F775W, F814W, F845M, FQ889N, FQ906N, F850LP, FQ924N, FQ937N, F953N</td>
        </tr>
        <tr>
          <th>ACS_HRC</th>
          <td>F220W, F250W, F330W, F344N, FR388N, F435W, FR459M, F475W, F502N, FR505N, F555W, F550M, F606W, F625W, FR656N, F658N, F660N, PR200L, POL_UV, POL_V, F775W, G800L, F814W, F892N, FR914M, F850LP</td>
        </tr>
        <tr>
          <th>HSP_POL</th>
          <td>F216M_0, F237M_0, F277M_0, F327M_0, F160LP_T, F160LP_A</td>
        </tr>
        <tr>
          <th>STIS_CCD</th>
          <td>F28X50LP_G230LB, 50CCD_G230LB, 50CORON_G230LB, F28X50LP_G230MB, 50CCD_G230MB, 50CORON_G230MB, F28X50OII_G430L, F28X50OII, F28X50OII_G430M, 50CCD_G430M, 50CORON_G430M, 50CCD_G430L, 50CORON_G430L, F28X50OIII_G430M, F28X50OIII, F28X50OIII_G430L, F28X50LP_G430M, F28X50LP_G430L, 50CCD, 50CORON, 50CCD_G750L, 50CORON_G750L, F28X50LP, F28X50LP_G750L, 50CCD_G750M, 50CORON_G750M, F28X50LP_G750M</td>
        </tr>
        <tr>
          <th>FOS_BLUE</th>
          <td>G130H, MIRROR, G190H, G160L, G400H, G270H, PRISM</td>
        </tr>
        <tr>
          <th>ACS_WFC</th>
          <td>FR388N, FR423N, F435W, FR459M, FR462N, F475W, F502N, FR505N, F555W, FR551N, F550M, FR601N, F606W, F625W, FR647M, FR656N, F658N, F660N, FR716N, POL_UV, POL_V, G800L, F775W, FR782N, F814W, FR853N, F892N, FR914M, F850LP, FR931N, FR1016N</td>
        </tr>
        <tr>
          <th>FOS_RED</th>
          <td>G780H, G190H, G160L, MIRROR, G270H, G400H, PRISM, G570H, G650L</td>
        </tr>
        <tr>
          <th>FGS</th>
          <td>F583W, ND5, PUPIL, F605W, F550W, F650W</td>
        </tr>
        <tr>
          <th>NICMOS1</th>
          <td>F090M, F095N, F097N, POL0S, POL240S, POL120S, F108N, F110M, F113N, F110W, F145M, F140W, F160W, F164N, F165M, F166N, F170M, F187N, F190N</td>
        </tr>
        <tr>
          <th>WFC3_IR</th>
          <td>F098M, G102, F105W, F110W, F125W, F126N, F127M, F128N, F130N, F132N, F139M, F140W, G141, F153M, F160W, F164N, F167N</td>
        </tr>
        <tr>
          <th>NICMOS3</th>
          <td>G096, F108N, F113N, F110W, G141, F150W, F160W, F164N, F166N, F187N, F175W, F190N, F196N, F200N, F212N, G206, F215N, F222M, F240M</td>
        </tr>
        <tr>
          <th>NICMOS2</th>
          <td>F110W, F160W, F165M, F171M, F180M, F187N, F187W, F190N, POL120L, POL0L, POL240L, F204M, F207M, F205W, F212N, F215N, F216N, F222M, F237M</td>
        </tr>
        <tr>
          <th>Integral</th>
          <th>OMC</th>
          <td>V_filter, V, V_opt</td>
        </tr>
        <tr>
          <th rowspan="3" valign="top">JWST</th>
          <th>NIRCam</th>
          <td>F070W, F090W, F115W, F140M, F150W, F162M, F164N, F150W2, F182M, F187N, F200W, F210M, F212N, F250M, F277W, F300M, F323N, F322W2, F335M, F356W, F360M, F405N, F410M, F430M, F444W, F460M, F466N, F470N, F480M</td>
        </tr>
        <tr>
          <th>NIRISS</th>
          <td>F090W, F115W, F140M, F150W, F158M, F200W, F277W, F356W, F380M, F430M, F444W, F480M</td>
        </tr>
        <tr>
          <th>MIRI</th>
          <td>F560W, F770W, F1000W, F1065C, F1140C, F1130W, F1280W, F1500W, F1550C, F1800W, F2100W, F2300C, F2550W</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">Keck</th>
          <th>LRIS</th>
          <td>NB4000, NB4325, B, g, NB5390, V, OG570, R, NB6741, Rs, GG495, I, NB8185, NB8560, RG850, NB9135, NB9148</td>
        </tr>
        <tr>
          <th>ESI</th>
          <td>B_fil, B, V, V_fil, R, R_fil</td>
        </tr>
        <tr>
          <th>NIRSPEC</th>
          <td>N1, N2, N3, N4, N5, FeII, thin, N6, H2, Kp, K, N7, KL, Lp, Mp, Mwide</td>
        </tr>
        <tr>
          <th>OSIRIS</th>
          <td>Zn2_spec, Zn3_spec, Zn3_imag, Zbb_spec, Zbb_imag, Zn4_spec, Zn5_spec, Jn1_spec, Jn1_imag, Jn2_imag, Jn2_spec, Jn3_imag, Jn3_spec, Jbb_spec, Jn4_spec, Hn1_spec, Hn1_imag, Hn2_imag, Hn2_spec, Hn3_spec, Hn3_imag, Hbb_spec, Hbb_imag, Hn4_spec, Hn4_imag, Hn5_imag, Hn5_spec, Kn1_imag, Kn1_spec, Kn2_spec, Kn2_imag, Kn3_spec, Kn3_imag, Kbb_spec, Kn4_spec, Kn4_imag, Kn5_imag, Kn5_spec</td>
        </tr>
        <tr>
          <th>NIRC2</th>
          <td>J, Hcont, H, Fell, Kp, Ks, Brgamma, K, Kcont, Lp, Ms</td>
        </tr>
        <tr>
          <th>LWS</th>
          <td>0800, 0890, 0990, Nwide, 1070, Spec10, 1170, SiC, 1250, 1765, 1875, 2310, 2450</td>
        </tr>
        <tr>
          <th>Kepler</th>
          <th>Kepler</th>
          <td>K</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">LBT</th>
          <th>LBCB</th>
          <td>bessel-U, bessel-B, bessel-V</td>
        </tr>
        <tr>
          <th>LBCR</th>
          <td>bessel-V, bessel-R, bessel-I</td>
        </tr>
        <tr>
          <th>LSST</th>
          <th>LSST</th>
          <td>u_filter, u, g_filter, g, r_filter, r, i, i_filter, z, z_filter, y, y_filter</td>
        </tr>
        <tr>
          <th>LaSilla</th>
          <th>SUSI2</th>
          <td>U, B_817, B, V, R, I, Z</td>
        </tr>
        <tr>
          <th>LasCumbres</th>
          <th>LasCumbres</th>
          <td>Bessel_B, Bessel_V, Bessel_R, Bessel_I</td>
        </tr>
        <tr>
          <th>McD</th>
          <th>DIAFI</th>
          <td>U, B, V, R, I</td>
        </tr>
        <tr>
          <th rowspan="3" valign="top">Misc</th>
          <th>MCPS</th>
          <td>U, B, V, I</td>
        </tr>
        <tr>
          <th>UCAC</th>
          <td>B, V</td>
        </tr>
        <tr>
          <th>APASS</th>
          <td>B, V</td>
        </tr>
        <tr>
          <th rowspan="3" valign="top">OSN</th>
          <th>Johnson</th>
          <td>U2, U, B2, B, V2, V, Cousins_R2, Cousins_R3, Cousins_R1, Cousins_I4, Cousins_I3, Cousins_I2, Cousins_I5, Cousins_I1</td>
        </tr>
        <tr>
          <th>Circ</th>
          <td>Johnson_B, Johnson_V, Cousins_R, Cousins_I</td>
        </tr>
        <tr>
          <th>Stromgren</th>
          <td>u2, u, v, v_z3, v_z2, v_z1, b, b_z1, b_z2, Hbetan, Hbetaw, Hbetan_z1, Hbetan_z2, Hbetaw_z1, Hbetaw_z2, Hbetaw_z3, y, y_z1, y_z2</td>
        </tr>
        <tr>
          <th>Paranal</th>
          <th>OmegaCAM</th>
          <td>B_aux_filter, B_qB_filter, B_qA_filter, B_qD_filter, B_filter, B_aux, B_qB, B_qA, B_qD, B, B_qC_filter, B_qC, V_qC, V_qB, V, V_qC_filter, V_qB_filter, V_qA, V_qD, V_filter, V_qA_filter, V_qD_filter, V_aux, V_aux_filter, v_strom_filter, v_strom, v_strom_aux_filter, v_strom_aux</td>
        </tr>
        <tr>
          <th>Roman</th>
          <th>WFI</th>
          <td>F062, F087, F106, F129, Prism, Grism, F146, F158, F184, F213</td>
        </tr>
        <tr>
          <th>SLOAN</th>
          <th>SDSS</th>
          <td>uprime_filter, u, g, gprime_filter, r, rprime_filter, i, iprime_filter, z, zprime_filter</td>
        </tr>
        <tr>
          <th>STELLA</th>
          <th>WiFSIP</th>
          <td>Strom_u, Strom_v, Strom_b, Strom_y</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">Scorpio</th>
          <th>Johnson</th>
          <td>U, B, V</td>
        </tr>
        <tr>
          <th>Cousins</th>
          <td>Rc, Ic</td>
        </tr>
        <tr>
          <th rowspan="3" valign="top">Spitzer</th>
          <th>IRAC</th>
          <td>I1, I2, I3, I4</td>
        </tr>
        <tr>
          <th>IRS</th>
          <td>Blue, Red</td>
        </tr>
        <tr>
          <th>MIPS</th>
          <td>24mu, 70mu, 160mu</td>
        </tr>
        <tr>
          <th>TESS</th>
          <th>TESS</th>
          <td>Red</td>
        </tr>
        <tr>
          <th>TJO</th>
          <th>MEIA</th>
          <td>U, B, V, Rc, Ic</td>
        </tr>
        <tr>
          <th>TNG</th>
          <th>TNG</th>
          <td>TNG01, TNG09, TNG10, TNG02, TNG11, TNG03, TNG04, TNG12, TNG13</td>
        </tr>
        <tr>
          <th>WFIRST</th>
          <th>WFI</th>
          <td>R062, Z087, Y106, J129, Prism, Grism, W146, H158, F184</td>
        </tr>
        <tr>
          <th>WISE</th>
          <th>WISE</th>
          <td>W1, W2, W3, W4</td>
        </tr>
      </tbody>
    </table>
