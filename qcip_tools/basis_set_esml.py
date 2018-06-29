"""
Basis sets obtained from ESML basis set exchange
"""

import requests
import bs4
from qcip_tools import atom as qcip_atom


class ESMLBasisSetError(Exception):
    pass


def base__get_atomic_basis_set(
        basis_set_name, path, atoms, js_peid, basis_set_format='Gaussian94', minimise=True):
    """Get basis set from ESML (low level version)

    :param basis_set_name: name of the basis set
    :type basis_set_name: str
    :param path: path in the ESML library
    :type path: str
    :param atoms: atoms requested (list of symbol, separated by a space)
    :type atoms: str
    :param basis_set_format: the format of the basis set
    :type basis_set_format: str
    :param minimise: optimized general contractions
    :type minimise: bool
    :param js_peid: the JS PEID
    :type js_peid: int
    :rtype: str
    """

    if basis_set_format not in AUTHORIZED_FORMATS:
        raise ESMLBasisSetError('type {} is unknown'.format(basis_set_format))

    if basis_set_name not in AVAILABLE_BS:
        raise ESMLBasisSetError('Basis {} is not available in the ESML'.format(basis_set_name))

    if not atoms:
        raise ESMLBasisSetError('empty list of atom')

    params = {
        'bsurl': path,
        'bsname': basis_set_name,
        'elts': atoms,
        'format': basis_set_format,
        'minimize': 'true' if minimise else 'false'
    }

    r = requests.get(URL.format(js_peid), params=params)
    if r.status_code != 200:
        raise ESMLBasisSetError('Status code is not 200')

    soup = bs4.BeautifulSoup(r.content, 'html.parser')
    r = soup.find('pre').string

    if 'EMSL  Basis Set Exchange Library' not in r:
        raise ESMLBasisSetError('Something wrong probably happen:\n {}'.format(r))

    return r


def get_atomic_basis_set(basis_set_name, atoms, basis_set_format='Gaussian94', minimise=True):
    """Get basis set from ESML (high level version)

    :param basis_set_name: name of the basis set
    :type basis_set_name: str
    :param atoms: atoms requested
    :type atoms: list
    :param basis_set_format: the format of the basis set
    :type basis_set_format: str
    :param minimise: optimized general contractions
    :type minimise: bool
    :rtype: str
    """

    if basis_set_name not in AVAILABLE_BS:
        raise ESMLBasisSetError('Basis {} is not available in the ESML'.format(basis_set_name))

    info = AVAILABLE_BS[basis_set_name]
    return info.get_basis(atoms, basis_set_format, minimise)


class ESMLBasisSet:
    """Store data for a given basis set

    :param name: name of the basis set in ESML
    :type name: str
    :param path: path to the basis set in the ESML
    :type path: str
    :param atoms: list of atoms for which the basis set is defined
    :type atoms: list
    :param typ: type of basis set (orbital, polarization, diffuse, ...)
    :type typ: str
    """

    def __init__(self, name, path, atoms, typ='orbital'):
        self.name = name
        self.path = path
        self.atoms = atoms
        self.type = typ

    def get_basis(self, atoms, basis_set_format='Gaussian94', minimise=True):
        """Get the basis set

        :param atoms: atoms requested
        :type atoms: list
        :param basis_set_format: the format of the basis set
        :type basis_set_format: str
        :param minimise: optimized general contractions
        :type minimise: bool
        :rtype: str
        """

        atoms_symbol = []
        for a in atoms:
            s = ''
            if type(a) is str:
                if a not in qcip_atom.SymbolToAtomicNumber:
                    raise ESMLBasisSetError('{} is not an atomic symbol'.format(a))
                s = a
            elif type(a) is int:
                if a not in qcip_atom.AtomicNumberToSymbol:
                    raise ESMLBasisSetError('{} is not a valid atomic number'.format(a))
                s = qcip_atom.AtomicNumberToSymbol[a]

            if s not in atoms_symbol:
                if s not in self.atoms:
                    raise ESMLBasisSetError('Basis set {} does not contains {}'.format(self.name, s))
                atoms_symbol.append(s)

        if not atoms_symbol:
            raise ESMLBasisSetError('empty list of atom')

        return base__get_atomic_basis_set(
            self.name,
            self.path,
            ' '.join(atoms_symbol),
            js_peid=JS_PEID,
            basis_set_format=basis_set_format,
            minimise=minimise)

# DO NOT EDIT BELOW THIS POINT!
# generated on December 18, 2017 (17:03)


#: current JS PEID
JS_PEID = 11543880926284

#: URL for the download
URL = 'https://bse.pnl.gov:443/bse/portal/user/anon/js_peid/{}/action/' \
      'portlets.BasisSetAction/template/courier_content/panel/Main/eventSubmit_doDownload/true'

#: List of authorized formats
AUTHORIZED_FORMATS = [
    'NWChem',
    'Gaussian94',
    'GAMESS-US',
    'GAMESS-UK',
    'Turbomole',
    'TX93',
    'Molpro',
    'MolproInt',
    'Hondo',
    'SuperMolecule',
    'Molcas',
    'HyperChem',
    'Dalton',
    'deMon-KS',
    'deMon2k',
    'AcesII'
]

# the basis sets:
AVAILABLE_BS = {
    'aug-cc-pwCVDZ-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pwCVDZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_39/AUG-CC-PWCVDZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'dhf-QZVP': ESMLBasisSet(
        'dhf-QZVP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/bert_link_5/DHF-QZVP.xml',
        ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pV(6+d)Z-RI': ESMLBasisSet(
        'cc-pV(6+d)Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_14/CC-PV6PDZ-RI.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Def2-SVP': ESMLBasisSet(
        'Def2-SVP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/weigend_new_0/DEF2-SVP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-pcJ-3_2006': ESMLBasisSet(
        'aug-pcJ-3_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_17/AUG-PCJ-3.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pV5Z_OPTRI': ESMLBasisSet(
        'aug-cc-pV5Z_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_22/AUG-CC-PV5Z_OPTRI.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVQZ-PP-RI': ESMLBasisSet(
        'cc-pVQZ-PP-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_20/CC-PVQZ-PP-RI.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'SARC2-QZVP-DKH/JK': ESMLBasisSet(
        'SARC2-QZVP-DKH/JK',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new_4/SARC2-QZVP-DKHJK.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
    'aug-cc-pV5Z': ESMLBasisSet(
        'aug-cc-pV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PV5Z-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pVQZ-DK': ESMLBasisSet(
        'aug-cc-pVQZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVQZ-DK-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'jul-cc-pV(T+d)Z': ESMLBasisSet(
        'jul-cc-pV(T+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_5/JUL-CC-PVTPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'SV + Double Rydberg (Dunning-Hay)': ESMLBasisSet(
        'SV + Double Rydberg (Dunning-Hay)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/SVPDOUBLERYDBERGDUNNING-HAY-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-cc-pVTZ-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pVTZ-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_27/AUG-CC-PVTZ-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'aug-cc-pwCV5Z-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pwCV5Z-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_42/AUG-CC-PWCV5Z-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'GAMESS PVTZ': ESMLBasisSet(
        'GAMESS PVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/GAMESSPVTZ-AGG.xml',
        ['H', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pcSseg-3': ESMLBasisSet(
        'aug-pcSseg-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_47/AUG-PCSSEG-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'DZQ': ESMLBasisSet(
        'DZQ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/amasunov_new/DZQ.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag']),  # noqa
    'cc-pVTZ(pt/sf/fw)': ESMLBasisSet(
        'cc-pVTZ(pt/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_PT_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-pcSseg-1': ESMLBasisSet(
        'aug-pcSseg-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_45/AUG-PCSSEG-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVQZ(pt/sf/lc)': ESMLBasisSet(
        'cc-pVQZ(pt/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_PT_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pwCVDZ-DK3': ESMLBasisSet(
        'cc-pwCVDZ-DK3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_36/CC-PWCVDZ-DK3.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']),  # noqa
    'maug-cc-pVDZ': ESMLBasisSet(
        'maug-cc-pVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Papajak_new_3/MAUG-CC-PVDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVQZ(fi/sf/fw)': ESMLBasisSet(
        'cc-pVQZ(fi/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_FI_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVTZ(seg-opt)': ESMLBasisSet(
        'cc-pVTZ(seg-opt)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_OPT1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'm6-31G': ESMLBasisSet(
        'm6-31G',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/marcel.swart_new/M6-31G.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']),  # noqa
    'Def2-QZVPPD': ESMLBasisSet(
        'Def2-QZVPPD',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/drappoport_new_3/DEF2-QZVPPD.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'Ahlrichs pVDZ': ESMLBasisSet(
        'Ahlrichs pVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AHLRICHSPVDZ-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pwCVQZ': ESMLBasisSet(
        'aug-cc-pwCVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCVQZ-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Br']),  # noqa
    'cc-pV6Z': ESMLBasisSet(
        'cc-pV6Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV6Z.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-31+G**': ESMLBasisSet(
        '6-31+G**',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_new_0/6-31PGSS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'IGLO-II': ESMLBasisSet(
        'IGLO-II',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/IGLOII.xml',
        ['H', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    's6-31G': ESMLBasisSet(
        's6-31G',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/marcel.swart_addElements_21/S6-31G.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'Pt - mDZP': ESMLBasisSet(
        'Pt - mDZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/diego_paschoal_new/PT-MDZP.xml',
        ['Pt']),  # noqa
    'cc-pCVDZ': ESMLBasisSet(
        'cc-pCVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCVDZ-AGG.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca']),  # noqa
    'cc-pCVQZ-F12': ESMLBasisSet(
        'cc-pCVQZ-F12',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_14/CC-PCVQZ-F12.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-311G(2df,2pd)': ESMLBasisSet(
        '6-311G(2df,2pd)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-311G2DF2PD-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'K', 'Ca']),  # noqa
    'aug-cc-pVTZ-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pVTZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_7/AUG-CC-PVTZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'NLO-V': ESMLBasisSet(
        'NLO-V',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/diego_paschoal_new_0/NLO-V.xml',
        ['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl']),  # noqa
    'Chipman DZP + Diffuse': ESMLBasisSet(
        'Chipman DZP + Diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CHIPMAN.xml',
        ['H', 'B', 'C', 'N', 'O', 'F']),  # noqa
    'cc-pVTZ(fi/sf/lc)': ESMLBasisSet(
        'cc-pVTZ(fi/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_FI_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    's3-21G*': ESMLBasisSet(
        's3-21G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/marcel.swart_addElements_12/S3-21GS.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'WTBS': ESMLBasisSet(
        'WTBS',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/WTBS.xml',
        ['He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'Weigend Coulomb Fitting': ESMLBasisSet(
        'Weigend Coulomb Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/edoapra_new_0/WEIGENDCOULOMBFITTING.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-cc-pV5Z-RI diffuse': ESMLBasisSet(
        'aug-cc-pV5Z-RI diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_10/AUG-CC-PV5Z-RIDIFFUSE.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'dhf-SV(P)': ESMLBasisSet(
        'dhf-SV(P)',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/bert_link_7/DHF-SVP.xml',
        ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    '3-21GSP': ESMLBasisSet(
        '3-21GSP',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/3-21GSP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pCV6Z (old)': ESMLBasisSet(
        'cc-pCV6Z (old)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCV6Z-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pCVQZ(old)': ESMLBasisSet(
        'cc-pCVQZ(old)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCVQZOLD-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcS-4': ESMLBasisSet(
        'pcS-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_28/PCS-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'NASA Ames cc-pCVQZ': ESMLBasisSet(
        'NASA Ames cc-pCVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_CC-PCVQZ.xml',
        ['Ti']),  # noqa
    'd-aug-cc-pV6Z': ESMLBasisSet(
        'd-aug-cc-pV6Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/D-AUG-CC-PV6Z-AGG.xml',
        ['H', 'B', 'C', 'N', 'O']),  # noqa
    'cc-pCV6Z': ESMLBasisSet(
        'cc-pCV6Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/edoapra_new_2/CC-PCV6Z_NEW.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pCVTZ-F12 MP2 Fitting': ESMLBasisSet(
        'cc-pCVTZ-F12 MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_53/CC-PCVTZ-F12MP2FITTING.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Lanl2-[5s4p4d2f]': ESMLBasisSet(
        'Lanl2-[5s4p4d2f]',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/omr_link_1/LANL2-[5S4P4D2F].xml',
        ['Rh']),  # noqa
    'cc-pCVDZ(old)': ESMLBasisSet(
        'cc-pCVDZ(old)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCVDZOLD-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pV5Z-DK': ESMLBasisSet(
        'aug-cc-pV5Z-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PV5Z-DK-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Binning/Curtiss VTZP': ESMLBasisSet(
        'Binning/Curtiss VTZP',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/BINNINGCURTISSVTZP-AGG.xml',
        ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Feller Misc. CVDZ': ESMLBasisSet(
        'Feller Misc. CVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/FELLER_VDZ.xml',
        ['K']),  # noqa
    'cc-pV5Z(pt/sf/fw)': ESMLBasisSet(
        'cc-pV5Z(pt/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z_PT_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '6-311G*': ESMLBasisSet(
        '6-311G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_addElements_2/6-311GS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'I']),  # noqa
    'cc-pwCVTZ-PP MP2 Fitting': ESMLBasisSet(
        'cc-pwCVTZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_23/CC-PWCVTZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'aug-mcc-pV5Z': ESMLBasisSet(
        'aug-mcc-pV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-MCC-PV5Z.xml',
        ['H']),  # noqa
    'maug-cc-pV(Q+d)Z': ESMLBasisSet(
        'maug-cc-pV(Q+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Papajak_new_1/MAUG-CC-PVQPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pcS-1': ESMLBasisSet(
        'aug-pcS-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_20/AUG-PCS-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pcSseg-0': ESMLBasisSet(
        'aug-pcSseg-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_44/AUG-PCSSEG-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'pcJ-0_2006': ESMLBasisSet(
        'pcJ-0_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_9/PCJ-0.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV(6+d)Z': ESMLBasisSet(
        'cc-pV(6+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV6Z_PLUS_D.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'SARC2-QZV-DKH/JK': ESMLBasisSet(
        'SARC2-QZV-DKH/JK',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new_7/SARC2-QZV-DKHJK.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
    'GAMESS VTZ': ESMLBasisSet(
        'GAMESS VTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/GAMESSVTZ.xml',
        ['H', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCVDZ-RI': ESMLBasisSet(
        'cc-pwCVDZ-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_28/CC-PWCVDZ-RI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    '5ZP-DKH': ESMLBasisSet(
        '5ZP-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_22/5ZP-DKH.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pwCV5Z': ESMLBasisSet(
        'aug-cc-pwCV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCV5Z-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVTZ-F12': ESMLBasisSet(
        'cc-pVTZ-F12',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_9/CC-PVTZ-F12.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVQZ-RI': ESMLBasisSet(
        'cc-pVQZ-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_4/CC-PVQZ-RI.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '3-21++G': ESMLBasisSet(
        '3-21++G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/3-21PPG-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'cc-pCV6Z(old)': ESMLBasisSet(
        'cc-pCV6Z(old)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCV6ZOLD-AGG.xml',
        ['O']),  # noqa
    'aug-cc-pVDZ-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pVDZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_6/AUG-CC-PVDZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'cc-pVQZ-PP MP2 Fitting': ESMLBasisSet(
        'cc-pVQZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_4/CC-PVQZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'MINI (Scaled)': ESMLBasisSet(
        'MINI (Scaled)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/MINISCALED.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'dhf-QZVPP': ESMLBasisSet(
        'dhf-QZVPP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/bert_link_6/DHF-QZVPP.xml',
        ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-cc-pVDZ-DK': ESMLBasisSet(
        'aug-cc-pVDZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVDZ-DK-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'x2c-TZVPPall-2c': ESMLBasisSet(
        'x2c-TZVPPall-2c',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_6/X2C-TZVPPALL-2C.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pwCVTZ-RI': ESMLBasisSet(
        'cc-pwCVTZ-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_29/CC-PWCVTZ-RI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'aug-cc-pVQZ': ESMLBasisSet(
        'aug-cc-pVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVQZ-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'd-aug-cc-pVDZ': ESMLBasisSet(
        'd-aug-cc-pVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/D-AUG-CC-PVDZ-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'ANO-RCC': ESMLBasisSet(
        'ANO-RCC',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/RolandLindh_new/ANO-RCC.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm']),  # noqa
    'NASA Ames cc-pV5Z': ESMLBasisSet(
        'NASA Ames cc-pV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_CC-PV5Z.xml',
        ['Ti', 'Fe']),  # noqa
    'maug-cc-pV(T+d)Z': ESMLBasisSet(
        'maug-cc-pV(T+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Papajak_new_0/MAUG-CC-PVTPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcJ-4_2006': ESMLBasisSet(
        'pcJ-4_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_13/PCJ-4.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Feller Misc. CVTZ': ESMLBasisSet(
        'Feller Misc. CVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/FELLER_VTZ.xml',
        ['K']),  # noqa
    'Lanl2-[10s8p7d3f2g]': ESMLBasisSet(
        'Lanl2-[10s8p7d3f2g]',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/omr_link_3/LANL2-[10S8P7D3F2G].xml',
        ['Rh']),  # noqa
    'aug-cc-pV6Z-RI diffuse': ESMLBasisSet(
        'aug-cc-pV6Z-RI diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_13/AUG-CC-PV6Z-RIDIFFUSE.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pVDZ_OPTRI': ESMLBasisSet(
        'aug-cc-pVDZ_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_19/AUG-CC-PVDZ_OPTRI.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCV5Z Core Set': ESMLBasisSet(
        'cc-pwCV5Z Core Set',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCV5Z0.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcS-0': ESMLBasisSet(
        'pcS-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_24/PCS-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'jul-cc-pV(D+d)Z': ESMLBasisSet(
        'jul-cc-pV(D+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new/JUL-CC-PVDPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pwCVTZ-NR MP2 Fitting': ESMLBasisSet(
        'aug-cc-pwCVTZ-NR MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_44/AUG-CC-PWCVTZ-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pV(D+d)Z': ESMLBasisSet(
        'cc-pV(D+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_PLUS_D.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcseg-3': ESMLBasisSet(
        'pcseg-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_32/PCSEG-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-pcJ-3': ESMLBasisSet(
        'aug-pcJ-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_13/AUG-PCJ-3.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcS-1': ESMLBasisSet(
        'pcS-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_25/PCS-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pcseg-3': ESMLBasisSet(
        'aug-pcseg-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_37/AUG-PCSEG-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-pcseg-4': ESMLBasisSet(
        'aug-pcseg-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_38/AUG-PCSEG-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVQZ(seg-opt)': ESMLBasisSet(
        'cc-pVQZ(seg-opt)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_OPT1.xml',
        ['H', 'He', 'C', 'O']),  # noqa
    'cc-pwCVTZ': ESMLBasisSet(
        'cc-pwCVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCVTZ-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pV7Z': ESMLBasisSet(
        'aug-pV7Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-PV7Z-AGG.xml',
        ['C', 'N', 'O', 'F', 'S']),  # noqa
    'aug-cc-pV5Z-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pV5Z-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_9/AUG-CC-PV5Z-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'Bauschlicher ANO': ESMLBasisSet(
        'Bauschlicher ANO',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/BAUSCLICHER_ANO.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']),  # noqa
    's3-21G': ESMLBasisSet(
        's3-21G',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/marcel.swart_addElements_3/S3-21G.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pCVTZ-F12_OPTRI': ESMLBasisSet(
        'cc-pCVTZ-F12_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_24/CC-PCVTZ-F12_OPTRI.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVDZ-F12_OPTRI': ESMLBasisSet(
        'cc-pVDZ-F12_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_16/CC-PVDZ-F12_OPTRI.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Wachters+f': ESMLBasisSet(
        'Wachters+f',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/WACHTERS.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']),  # noqa
    'may-cc-pV(T+d)Z': ESMLBasisSet(
        'may-cc-pV(T+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_18/MAY-CC-PVTPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Def2-QZVP': ESMLBasisSet(
        'Def2-QZVP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/weigend_new_0/DEF2-QZVP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pCV5Z': ESMLBasisSet(
        'cc-pCV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_7/CC-PCV5Z.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca']),  # noqa
    'pc-1': ESMLBasisSet(
        'pc-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_1/PC-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pCVQZ-F12_OPTRI': ESMLBasisSet(
        'cc-pCVQZ-F12_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_25/CC-PCVQZ-F12_OPTRI.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVTZ-DK': ESMLBasisSet(
        'cc-pVTZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ-DK.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'pcJ-3': ESMLBasisSet(
        'pcJ-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_8/PCJ-3.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCV5Z': ESMLBasisSet(
        'cc-pwCV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCV5Z-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pCVQZ': ESMLBasisSet(
        'cc-pCVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCVQZ-AGG.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca']),  # noqa
    'pcS-2': ESMLBasisSet(
        'pcS-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_26/PCS-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcSseg-4': ESMLBasisSet(
        'pcSseg-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_43/PCSSEG-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '6-31G*-Blaudeau': ESMLBasisSet(
        '6-31G*-Blaudeau',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-31GS-BLAUDEAU-AGG.xml',
        ['K', 'Ca']),  # noqa
    'aug-cc-pwCV5Z-NR': ESMLBasisSet(
        'aug-cc-pwCV5Z-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCV5Z-NR-AGG.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'B2 basis set for Zn': ESMLBasisSet(
        'B2 basis set for Zn',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/eamin_new/MODIFIEDWACHTERS-HAYBASISFORZN.xml',
        ['Zn']),  # noqa
    '6-311++G*': ESMLBasisSet(
        '6-311++G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_new_4/6-311PPGS.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    '6ZP': ESMLBasisSet(
        '6ZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_10/6ZP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'dhf-TZVPP': ESMLBasisSet(
        'dhf-TZVPP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/bert_link_10/DHF-TZVPP.xml',
        ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-cc-pV(Q+d)Z': ESMLBasisSet(
        'aug-cc-pV(Q+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVQPDZ-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pVTZ': ESMLBasisSet(
        'aug-cc-pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVTZ-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVTZ-F12 MP2 Fitting': ESMLBasisSet(
        'cc-pVTZ-F12 MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_50/CC-PVTZ-F12MP2FITTING.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6ZaPa-NR': ESMLBasisSet(
        '6ZaPa-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_1/6ZAPA-NR.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '2ZaPa-NR': ESMLBasisSet(
        '2ZaPa-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_6/2ZAPA-NR.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'ATZP': ESMLBasisSet(
        'ATZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_15/ATZP.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']),  # noqa
    'QZP': ESMLBasisSet(
        'QZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_13/QZP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']),  # noqa
    'pc-3': ESMLBasisSet(
        'pc-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_3/PC-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Lanl2-[6s4p4d2f]': ESMLBasisSet(
        'Lanl2-[6s4p4d2f]',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/omr_link_2/LANL2-[6S4P4D2F].xml',
        ['Rh']),  # noqa
    'aug-cc-pVDZ-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pVDZ-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_26/AUG-CC-PVDZ-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'aug-mcc-pV7Z': ESMLBasisSet(
        'aug-mcc-pV7Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-MCC-PV7Z.xml',
        ['H']),  # noqa
    'aug-pcJ-2': ESMLBasisSet(
        'aug-pcJ-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_12/AUG-PCJ-2.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pwCVTZ-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pwCVTZ-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_31/AUG-CC-PWCVTZ-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'Ahlrichs TZV': ESMLBasisSet(
        'Ahlrichs TZV',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AHLRICHS_TZV.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-pcJ-4_2006': ESMLBasisSet(
        'aug-pcJ-4_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_18/AUG-PCJ-4.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'x2c-coulomb-fitting': ESMLBasisSet(
        'x2c-coulomb-fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_7/X2C-COULOMB-FITTING.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'DZP-DKH': ESMLBasisSet(
        'DZP-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_19/DZP-DKH.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-cc-pVQZ-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pVQZ-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_28/AUG-CC-PVQZ-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'aug-cc-pV(6+d)Z': ESMLBasisSet(
        'aug-cc-pV(6+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PV6PDZ-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-311+G*-J': ESMLBasisSet(
        '6-311+G*-J',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/sauer_new_2/6-311PGS-J.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    'cc-pVTZ-F12_OPTRI': ESMLBasisSet(
        'cc-pVTZ-F12_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_17/CC-PVTZ-F12_OPTRI.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV(T+d)Z-DK': ESMLBasisSet(
        'cc-pV(T+d)Z-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/wanyi_new_0/CC-PVTPDZ-DK.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'coemd-4': ESMLBasisSet(
        'coemd-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_0/COEMD-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'DZP': ESMLBasisSet(
        'DZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_17/DZP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'LANL2TZ': ESMLBasisSet(
        'LANL2TZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/linnyr419_link/LANL2TZ.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']),  # noqa
    'pcSseg-0': ESMLBasisSet(
        'pcSseg-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_39/PCSSEG-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pVDZ-PP-RI Diffuse': ESMLBasisSet(
        'aug-cc-pVDZ-PP-RI Diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_22/AUG-CC-PVDZ-PP-RIDIFFUSE.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    '4ZaPa-NR_CV': ESMLBasisSet(
        '4ZaPa-NR_CV',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_17/4ZAPA-NR_CV.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVDZ-DK3': ESMLBasisSet(
        'cc-pVDZ-DK3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_30/CC-PVDZ-DK3.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']),  # noqa
    'NASA Ames ANO': ESMLBasisSet(
        'NASA Ames ANO',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_ANO.xml',
        ['H', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'P', 'Ti', 'Fe', 'Ni']),  # noqa
    'cc-pV(Q+d)Z-RI': ESMLBasisSet(
        'cc-pV(Q+d)Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_2/CC-PVQPDZ-RI.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '3-21++G*': ESMLBasisSet(
        '3-21++G*',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/3-21PPGS-AGG.xml',
        ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVTZ MP2 Fitting': ESMLBasisSet(
        'cc-pVTZ MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_1/CC-PVTZ-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    '3-21G*': ESMLBasisSet(
        '3-21G*',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/3-21GS-AGG.xml',
        ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-311G-J': ESMLBasisSet(
        '6-311G-J',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/sauer_correction_2/6-311G-J.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    'aug-cc-pV(7+d)Z': ESMLBasisSet(
        'aug-cc-pV(7+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_22/AUG-CC-PV7PDZ.xml',
        ['Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    'aug-cc-pwCVQZ-DK': ESMLBasisSet(
        'aug-cc-pwCVQZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCVQZ-DK-AGG.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'aug-cc-pV5Z-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pV5Z-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_29/AUG-CC-PV5Z-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'aug-pcJ-2_2006': ESMLBasisSet(
        'aug-pcJ-2_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_16/AUG-PCJ-2.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pV7Z': ESMLBasisSet(
        'aug-cc-pV7Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_23/AUG-CC-PV7Z.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'cc-pwCVDZ-RI tight': ESMLBasisSet(
        'cc-pwCVDZ-RI tight',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_15/CC-PWCVDZ-RITIGHT.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'A5ZP': ESMLBasisSet(
        'A5ZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_12/A5ZP.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pVQZ-PP-RI Diffuse': ESMLBasisSet(
        'aug-cc-pVQZ-PP-RI Diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_24/AUG-CC-PVQZ-PP-RIDIFFUSE.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'x2c-TZVPall': ESMLBasisSet(
        'x2c-TZVPall',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_1/X2C-TZVPALL.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'NMR-DKH (TZ2P)': ESMLBasisSet(
        'NMR-DKH (TZ2P)',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/diego_paschoal_new_1/NMR-DKHTZ2P.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Pt']),  # noqa
    'pcJ-4': ESMLBasisSet(
        'pcJ-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_9/PCJ-4.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pVTZ-DK': ESMLBasisSet(
        'aug-cc-pVTZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVTZ-DK-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'cc-pV5Z(pt/sf/sc)': ESMLBasisSet(
        'cc-pV5Z(pt/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z_PT_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pwCVTZ-DK': ESMLBasisSet(
        'cc-pwCVTZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCVTZ-DK.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'cc-pVQZ-DK': ESMLBasisSet(
        'cc-pVQZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ-DK.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '6-311+G**': ESMLBasisSet(
        '6-311+G**',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_new_2/6-311PGSS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'aug-cc-pVTZ-PP-RI Diffuse': ESMLBasisSet(
        'aug-cc-pVTZ-PP-RI Diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_23/AUG-CC-PVTZ-PP-RIDIFFUSE.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    '3ZaPa-NR': ESMLBasisSet(
        '3ZaPa-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_5/3ZAPA-NR.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCVTZ-RI tight': ESMLBasisSet(
        'cc-pwCVTZ-RI tight',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_16/CC-PWCVTZ-RITIGHT.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Cologne DKH2': ESMLBasisSet(
        'Cologne DKH2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/mdolg_new/COLOGNEDKH2.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']),  # noqa
    'MINIs-BSIP1': ESMLBasisSet(
        'MINIs-BSIP1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/aoterodelaroza_correction/MINIS-BSIP1.xml',
        ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl']),  # noqa
    'SBKJC Polarized (p,2d) - LFK': ESMLBasisSet(
        'SBKJC Polarized (p,2d) - LFK',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/nick.labello_link/SBKJCPOLARIZEDP2D-LFK.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'pcemd-3': ESMLBasisSet(
        'pcemd-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_4/PCEMD-3.xml',
        ['He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar', 'P', 'Cl', 'Al', 'Si', 'S']),  # noqa
    'MG3S': ESMLBasisSet(
        'MG3S',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/lever046_new/MG3S.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'STO-2G': ESMLBasisSet(
        'STO-2G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/STO-2G.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sr']),  # noqa
    'cc-pCVDZ-F12 MP2 Fitting': ESMLBasisSet(
        'cc-pCVDZ-F12 MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_52/CC-PCVDZ-F12MP2FITTING.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pVTZ_OPTRI': ESMLBasisSet(
        'aug-cc-pVTZ_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_20/AUG-CC-PVTZ_OPTRI.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCVTZ-NR': ESMLBasisSet(
        'cc-pwCVTZ-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCVTZ-NR.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'aug-pc-1': ESMLBasisSet(
        'aug-pc-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_7/AUG-PC-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'AQZP': ESMLBasisSet(
        'AQZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_14/AQZP.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']),  # noqa
    'm6-31G*': ESMLBasisSet(
        'm6-31G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/marcel.swart_new_14/M6-31GS.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']),  # noqa
    'SV (Dunning-Hay)': ESMLBasisSet(
        'SV (Dunning-Hay)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DUNHAYSV.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'cc-pV(T+d)Z-RI': ESMLBasisSet(
        'cc-pV(T+d)Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_1/CC-PVTPDZ-RI.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '5ZP': ESMLBasisSet(
        '5ZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_11/5ZP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCVTZ-NR MP2 Fitting': ESMLBasisSet(
        'cc-pwCVTZ-NR MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_43/CC-PWCVTZ-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pVDZ(seg-opt)': ESMLBasisSet(
        'cc-pVDZ(seg-opt)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_OPT1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'pV6Z': ESMLBasisSet(
        'pV6Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV6Z-OLD.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    'cc-pVQZ(pt/sf/fw)': ESMLBasisSet(
        'cc-pVQZ(pt/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_PT_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Binning/Curtiss VTZ': ESMLBasisSet(
        'Binning/Curtiss VTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/BINCURT-VTZ.xml',
        ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pwCVQZ-NR': ESMLBasisSet(
        'cc-pwCVQZ-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCVQZ-NR.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'SVP (Dunning-Hay)': ESMLBasisSet(
        'SVP (Dunning-Hay)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/SVPDUNNING-HAY-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'pcJ-2': ESMLBasisSet(
        'pcJ-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_7/PCJ-2.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-31G**': ESMLBasisSet(
        '6-31G**',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-31GSS-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'NASA Ames cc-pCV5Z': ESMLBasisSet(
        'NASA Ames cc-pCV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_CC-PCV5Z.xml',
        ['Ti']),  # noqa
    'cc-pCVTZ(old)': ESMLBasisSet(
        'cc-pCVTZ(old)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCVTZOLD-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV(D+d)Z-RI': ESMLBasisSet(
        'cc-pV(D+d)Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_0/CC-PVDPDZ-RI.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVDZ-RI': ESMLBasisSet(
        'cc-pVDZ-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new/CC-PVDZ-RI.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pwCVQZ-PP MP2 Fitting': ESMLBasisSet(
        'cc-pwCVQZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_24/CC-PWCVQZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'aug-cc-pVTZ-J': ESMLBasisSet(
        'aug-cc-pVTZ-J',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/sauer_correction_0/AUG-CC-PVTZ-J.xml',
        ['H', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Se']),  # noqa
    'cc-pV6Z-RI': ESMLBasisSet(
        'cc-pV6Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_12/CC-PV6Z-RI.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Def2-TZVP': ESMLBasisSet(
        'Def2-TZVP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/weigend_new_0/DEF2-TZVP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'Def2-SV(P)': ESMLBasisSet(
        'Def2-SV(P)',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/weigend_new_0/DEF2-SVPP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'Feller Misc. CVQZ': ESMLBasisSet(
        'Feller Misc. CVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/FELLER_VQZ.xml',
        ['K']),  # noqa
    'TZP-DKH': ESMLBasisSet(
        'TZP-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_20/TZP-DKH.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']),  # noqa
    'pcSseg-1': ESMLBasisSet(
        'pcSseg-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_40/PCSSEG-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pwCV5Z-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pwCV5Z-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_33/AUG-CC-PWCV5Z-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'SV + Rydberg (Dunning-Hay)': ESMLBasisSet(
        'SV + Rydberg (Dunning-Hay)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/SVPRYDBERGDUNNING-HAY-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-pcJ-0': ESMLBasisSet(
        'aug-pcJ-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_10/AUG-PCJ-0.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV9Z': ESMLBasisSet(
        'cc-pV9Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV9Z.xml',
        ['Ne']),  # noqa
    'cc-pwCVQZ': ESMLBasisSet(
        'cc-pwCVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCVQZ-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Br']),  # noqa
    's6-31G*': ESMLBasisSet(
        's6-31G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/marcel.swart_addElements_33/S6-31GS.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'aug-pcS-2': ESMLBasisSet(
        'aug-pcS-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_21/AUG-PCS-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'ccJ-pV5Z': ESMLBasisSet(
        'ccJ-pV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/beud_new_2/CCJ-PV5Z.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'cc-pVDZ': ESMLBasisSet(
        'cc-pVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'pcJ-2_2006': ESMLBasisSet(
        'pcJ-2_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_11/PCJ-2.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '3-21G': ESMLBasisSet(
        '3-21G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/3-21G.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs']),  # noqa
    'cc-pVDZ(pt/sf/sc)': ESMLBasisSet(
        'cc-pVDZ(pt/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_PT_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'NASA Ames cc-pVQZ': ESMLBasisSet(
        'NASA Ames cc-pVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_CC-PVQZ.xml',
        ['Ti', 'Fe']),  # noqa
    'DZP + Rydberg (Dunning)': ESMLBasisSet(
        'DZP + Rydberg (Dunning)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DZPPRYDBERGDUNNING-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    '6-31G(2df,p)': ESMLBasisSet(
        '6-31G(2df,p)',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/edoapra_new_1/6-31G2DFP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pwCVQZ-NR': ESMLBasisSet(
        'aug-cc-pwCVQZ-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCVQZ-NR-AGG.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pwCV5Z-DK': ESMLBasisSet(
        'cc-pwCV5Z-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCV5Z-DK.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'Binning/Curtiss SVP': ESMLBasisSet(
        'Binning/Curtiss SVP',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/BINNINGCURTISSSVP-AGG.xml',
        ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'MIDI (Huzinaga)': ESMLBasisSet(
        'MIDI (Huzinaga)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/MIDI.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Cs']),  # noqa
    'Def2-TZVPPD': ESMLBasisSet(
        'Def2-TZVPPD',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/drappoport_new_1/DEF2-TZVPPD.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pVQZ-F12_OPTRI': ESMLBasisSet(
        'cc-pVQZ-F12_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_18/CC-PVQZ-F12_OPTRI.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Def2-QZVPP': ESMLBasisSet(
        'Def2-QZVPP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/weigend_new_0/DEF2-QZVPP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pVQZ(pt/sf/sc)': ESMLBasisSet(
        'cc-pVQZ(pt/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_PT_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-pc-3': ESMLBasisSet(
        'aug-pc-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction/AUG-PC-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pVDZ-RI diffuse': ESMLBasisSet(
        'aug-cc-pVDZ-RI diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_7/AUG-CC-PVDZ-RIDIFFUSE.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'pc-4': ESMLBasisSet(
        'pc-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_4/PC-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pVQZ-NR MP2 Fitting': ESMLBasisSet(
        'aug-cc-pVQZ-NR MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_47/AUG-CC-PVQZ-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'SARC2-QZV-ZORA/JK': ESMLBasisSet(
        'SARC2-QZV-ZORA/JK',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new_6/SARC2-QZV-ZORAJK.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
    'Partridge Uncontracted 3': ESMLBasisSet(
        'Partridge Uncontracted 3',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/PARTRIDGE3.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'pcseg-0': ESMLBasisSet(
        'pcseg-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_29/PCSEG-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'LANL2TZ+': ESMLBasisSet(
        'LANL2TZ+',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/linnyr419_link_4/LANL2TZP.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pVTZ-PP-RI': ESMLBasisSet(
        'cc-pVTZ-PP-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_19/CC-PVTZ-PP-RI.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pwCVDZ-PP MP2 Fitting': ESMLBasisSet(
        'cc-pwCVDZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_22/CC-PWCVDZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'pc-0': ESMLBasisSet(
        'pc-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_0/PC-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '6-311+G': ESMLBasisSet(
        '6-311+G',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_new_1/6-311PG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'Lanl2DZ+1d1f': ESMLBasisSet(
        'Lanl2DZ+1d1f',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/omr_link/LANL2DZP1D1F.xml',
        ['Rh']),  # noqa
    'Roos Augmented Triple Zeta ANO': ESMLBasisSet(
        'Roos Augmented Triple Zeta ANO',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/ROOS_ANO2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'aug-pcseg-0': ESMLBasisSet(
        'aug-pcseg-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_34/AUG-PCSEG-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pCVTZ': ESMLBasisSet(
        'aug-cc-pCVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PCVTZ-AGG.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVQZ-DK3': ESMLBasisSet(
        'cc-pVQZ-DK3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_35/CC-PVQZ-DK3.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']),  # noqa
    '5ZaPa-NR': ESMLBasisSet(
        '5ZaPa-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_2/5ZAPA-NR.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'd-aug-cc-pVTZ': ESMLBasisSet(
        'd-aug-cc-pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/D-AUG-CC-PVTZ-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-cc-pCV(T+d)Z': ESMLBasisSet(
        'aug-cc-pCV(T+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/wanyi_new/AUG-CC-PCVTPDZ.xml',
        ['Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    'pcseg-4': ESMLBasisSet(
        'pcseg-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_33/PCSEG-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-pcS-4': ESMLBasisSet(
        'aug-pcS-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_23/AUG-PCS-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcemd-2': ESMLBasisSet(
        'pcemd-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_3/PCEMD-2.xml',
        ['He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar', 'P', 'Cl', 'Al', 'Si', 'S']),  # noqa
    'Roos Augmented Double Zeta ANO': ESMLBasisSet(
        'Roos Augmented Double Zeta ANO',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/ROOS_ANO1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    '6-31G': ESMLBasisSet(
        '6-31G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-31G.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'aug-cc-pVQZ_OPTRI': ESMLBasisSet(
        'aug-cc-pVQZ_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_21/AUG-CC-PVQZ_OPTRI.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-31++G**': ESMLBasisSet(
        '6-31++G**',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-31PPGSS-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    '6-31++G**-J': ESMLBasisSet(
        '6-31++G**-J',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/sauer_new_5/6-31PPGSS-J.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    '6-31++G': ESMLBasisSet(
        '6-31++G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-31PPG-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'cc-pVDZ-F12': ESMLBasisSet(
        'cc-pVDZ-F12',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_8/CC-PVDZ-F12.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Def2-QZVPD': ESMLBasisSet(
        'Def2-QZVPD',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/drappoport_new_2/DEF2-QZVPD.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'SVP + Diffuse (Dunning-Hay)': ESMLBasisSet(
        'SVP + Diffuse (Dunning-Hay)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/SVPPDIFFUSEDUNNING-HAY-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'DZP + Diffuse (Dunning)': ESMLBasisSet(
        'DZP + Diffuse (Dunning)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DZPPDIFFUSEDUNNING-AGG.xml',
        ['H', 'Li', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-pcJ-0_2006': ESMLBasisSet(
        'aug-pcJ-0_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_14/AUG-PCJ-0.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pwCVTZ-NR': ESMLBasisSet(
        'aug-cc-pwCVTZ-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCVTZ-NR-AGG.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'dhf-SVP': ESMLBasisSet(
        'dhf-SVP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/bert_link_8/DHF-SVP.xml',
        ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pVTZ(fi/sf/fw)': ESMLBasisSet(
        'cc-pVTZ(fi/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_FI_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pV5Z-RI': ESMLBasisSet(
        'cc-pV5Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_5/CC-PV5Z-RI.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'maug-cc-pV(D+d)Z': ESMLBasisSet(
        'maug-cc-pV(D+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Papajak_new_2/MAUG-CC-PVDPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'x2c-SVPall': ESMLBasisSet(
        'x2c-SVPall',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_0/X2C-SVPALL.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    '6-31G-J': ESMLBasisSet(
        '6-31G-J',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/sauer_correction_1/6-31G-J.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    'McLean/Chandler VTZ': ESMLBasisSet(
        'McLean/Chandler VTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/MCLEAN-CHAND.xml',
        ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-31G*': ESMLBasisSet(
        '6-31G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_addElements_1/6-31GS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'IGLO-III': ESMLBasisSet(
        'IGLO-III',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/IGLOIII.xml',
        ['H', 'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    'aug-cc-pwCVTZ-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pwCVTZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_40/AUG-CC-PWCVTZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'cc-pVQZ-NR MP2 Fitting': ESMLBasisSet(
        'cc-pVQZ-NR MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_45/CC-PVQZ-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'un-pcemd-ref': ESMLBasisSet(
        'un-pcemd-ref',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_9/UN-PCEMD-REF.xml',
        ['He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar', 'P', 'Cl', 'Al', 'Si', 'S']),  # noqa
    'aug-pcS-3': ESMLBasisSet(
        'aug-pcS-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_22/AUG-PCS-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCVQZ-DK3': ESMLBasisSet(
        'cc-pwCVQZ-DK3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_38/CC-PWCVQZ-DK3.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']),  # noqa
    '6-31++G*': ESMLBasisSet(
        '6-31++G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_addElements/6-31PPGS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'cc-pVQZ(fi/sf/lc)': ESMLBasisSet(
        'cc-pVQZ(fi/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_FI_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '6-311G**': ESMLBasisSet(
        '6-311G**',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-311GSS-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'I']),  # noqa
    'ccJ-pVQZ': ESMLBasisSet(
        'ccJ-pVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/beud_new_1/CCJ-PVQZ.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-cc-pV5Z-PP-RI Diffuse': ESMLBasisSet(
        'aug-cc-pV5Z-PP-RI Diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_26/AUG-CC-PV5Z-PP-RIDIFFUSE.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pV5Z(fi/sf/sc)': ESMLBasisSet(
        'cc-pV5Z(fi/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z_FI_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'coemd-3': ESMLBasisSet(
        'coemd-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_1/COEMD-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcemd-4': ESMLBasisSet(
        'pcemd-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_5/PCEMD-4.xml',
        ['He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar', 'P', 'Cl', 'Al', 'Si', 'S']),  # noqa
    'DZ + Rydberg (Dunning)': ESMLBasisSet(
        'DZ + Rydberg (Dunning)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DZPRYDBERGDUNNING-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    '6-311++G(2d,2p)': ESMLBasisSet(
        '6-311++G(2d,2p)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-311PPG2D2P-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'Def2-TZVPP': ESMLBasisSet(
        'Def2-TZVPP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/weigend_new_0/DEF2-TZVPP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pwCVQZ-RI tight': ESMLBasisSet(
        'cc-pwCVQZ-RI tight',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_17/CC-PWCVQZ-RITIGHT.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pV(5+d)Z': ESMLBasisSet(
        'aug-cc-pV(5+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PV5PDZ-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'LANL08d': ESMLBasisSet(
        'LANL08d',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/linnyr419_link_5/LANL08D.xml',
        ['Si', 'P', 'S', 'Cl', 'Ge', 'As', 'Se', 'Br', 'Sn', 'Sb', 'Te', 'I', 'Pb', 'Bi']),  # noqa
    'd-aug-cc-pV5Z': ESMLBasisSet(
        'd-aug-cc-pV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/D-AUG-CC-PV5Z-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    '5ZaPa-NR_CV': ESMLBasisSet(
        '5ZaPa-NR_CV',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_20/5ZAPA-NR_CV.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-311+G*': ESMLBasisSet(
        '6-311+G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_addElements_3/6-311PGS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'cc-pCVDZ-F12_OPTRI': ESMLBasisSet(
        'cc-pCVDZ-F12_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_23/CC-PCVDZ-F12_OPTRI.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'CVTZ': ESMLBasisSet(
        'CVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/ErinDahlke_new/CVTZ_LIBENAMGKCA.xml',
        ['Li', 'Be', 'Na', 'Mg', 'K', 'Ca']),  # noqa
    '6-31+G*': ESMLBasisSet(
        '6-31+G*',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_addElements_0/6-31PGS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'cc-pCVTZ': ESMLBasisSet(
        'cc-pCVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PCVTZ-AGG.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca']),  # noqa
    'LANL08': ESMLBasisSet(
        'LANL08',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/linnyr419_link_0/LANL08.xml',
        ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']),  # noqa
    '6-31+G': ESMLBasisSet(
        '6-31+G',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_new/6-31PG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'cc-pV(Q+d)Z': ESMLBasisSet(
        'cc-pV(Q+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_PLUS_D.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'SVP + Rydberg (Dunning-Hay)': ESMLBasisSet(
        'SVP + Rydberg (Dunning-Hay)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/SVPPRYDBERGDUNNING-HAY-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'coemd-ref': ESMLBasisSet(
        'coemd-ref',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new/COEMD-REF.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Sadlej+': ESMLBasisSet(
        'Sadlej+',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/DFT_new/SADLEJP.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    '3ZaPa-NR_CV': ESMLBasisSet(
        '3ZaPa-NR_CV',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_10/3ZAPA-NR_CV.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcseg-2': ESMLBasisSet(
        'pcseg-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_31/PCSEG-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVTZ-DK3': ESMLBasisSet(
        'cc-pVTZ-DK3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_34/CC-PVTZ-DK3.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']),  # noqa
    '7ZaPa-NR': ESMLBasisSet(
        '7ZaPa-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_7/7ZAPA-NR.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'jun-cc-pV(T+d)Z': ESMLBasisSet(
        'jun-cc-pV(T+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_17/JUN-CC-PVTPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-mcc-pV8Z': ESMLBasisSet(
        'aug-mcc-pV8Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-MCC-PV8Z.xml',
        ['H']),  # noqa
    'cc-pwCVDZ': ESMLBasisSet(
        'cc-pwCVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCVDZ-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV5Z-PP-RI': ESMLBasisSet(
        'cc-pV5Z-PP-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_27/CC-PV5Z-PP-RI.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'cc-pVQZ-F12': ESMLBasisSet(
        'cc-pVQZ-F12',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_10/CC-PVQZ-F12.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'DZP (Dunning)': ESMLBasisSet(
        'DZP (Dunning)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DZPDUNNING-AGG.xml',
        ['H', 'Li', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    'x2c-TZVPall-2c': ESMLBasisSet(
        'x2c-TZVPall-2c',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_5/X2C-TZVPALL-2C.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'modified LANL2DZ': ESMLBasisSet(
        'modified LANL2DZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/cewebstr_link/MODIFIEDLANL2DZ.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']),  # noqa
    'dhf-TZVP': ESMLBasisSet(
        'dhf-TZVP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/bert_link_9/DHF-TZVP.xml',
        ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'x2c-SV(P)all': ESMLBasisSet(
        'x2c-SV(P)all',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new/X2C-SVPALL.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-cc-pV5Z-NR MP2 Fitting': ESMLBasisSet(
        'aug-cc-pV5Z-NR MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_48/AUG-CC-PV5Z-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'x2c-TZVPPall': ESMLBasisSet(
        'x2c-TZVPPall',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_2/X2C-TZVPPALL.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    '6-31G(3df,3pd)': ESMLBasisSet(
        '6-31G(3df,3pd)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-31G3DF3PD-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV(5+d)Z': ESMLBasisSet(
        'cc-pV(5+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z_PLUS_D.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'un-ccemd-ref': ESMLBasisSet(
        'un-ccemd-ref',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_8/UN-CCEMD-REF.xml',
        ['He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar', 'P', 'Cl', 'Al', 'Si', 'S']),  # noqa
    'SVP + Diffuse + Rydberg': ESMLBasisSet(
        'SVP + Diffuse + Rydberg',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/SVPPDIFFUSEPRYDBERG-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'LANL08+': ESMLBasisSet(
        'LANL08+',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/linnyr419_link_3/LANL08P.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pVDZ(pt/sf/fw)': ESMLBasisSet(
        'cc-pVDZ(pt/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_PT_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pV(5+d)Z-RI': ESMLBasisSet(
        'cc-pV(5+d)Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_6/CC-PV5PDZ-RI.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'UGBS': ESMLBasisSet(
        'UGBS',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/edoapra_new/UGBS.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pu', 'Am', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']),  # noqa
    'Sadlej pVTZ': ESMLBasisSet(
        'Sadlej pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/SADLEJ.xml',
        ['H', 'Li', 'Be', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Si', 'P', 'S', 'Cl', 'K', 'Ca', 'Br', 'Rb', 'Sr', 'I']),  # noqa
    'ccJ-pVTZ': ESMLBasisSet(
        'ccJ-pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/beud_new_0/CCJ-PVTZ.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-pcseg-1': ESMLBasisSet(
        'aug-pcseg-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_35/AUG-PCSEG-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'maug-cc-pVQZ': ESMLBasisSet(
        'maug-cc-pVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Papajak_new_6/MAUG-CC-PVQZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV5Z-PP MP2 Fitting': ESMLBasisSet(
        'cc-pV5Z-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_5/CC-PV5Z-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'aug-cc-pCVDZ': ESMLBasisSet(
        'aug-cc-pCVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PCVDZ-AGG.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-311++G(3df,3pd)': ESMLBasisSet(
        '6-311++G(3df,3pd)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-311PPG3DF3PD-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcJ-1': ESMLBasisSet(
        'pcJ-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_6/PCJ-1.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-31G-Blaudeau': ESMLBasisSet(
        '6-31G-Blaudeau',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-31G-BLAUDEAU.xml',
        ['K', 'Ca']),  # noqa
    'd-aug-cc-pVQZ': ESMLBasisSet(
        'd-aug-cc-pVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/D-AUG-CC-PVQZ-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-cc-pwCVDZ-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pwCVDZ-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_30/AUG-CC-PWCVDZ-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'cc-pVTZ(pt/sf/sc)': ESMLBasisSet(
        'cc-pVTZ(pt/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_PT_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVTZ': ESMLBasisSet(
        'cc-pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Def2-TZVPD': ESMLBasisSet(
        'Def2-TZVPD',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/drappoport_new_0/DEF2-TZVPD.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'G3MP2LargeXP': ESMLBasisSet(
        'G3MP2LargeXP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_25/G3MP2LARGEXP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVDZ(pt/sf/lc)': ESMLBasisSet(
        'cc-pVDZ(pt/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_PT_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVDZ-DK': ESMLBasisSet(
        'cc-pVDZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ-DK.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'LANL2TZ(f)': ESMLBasisSet(
        'LANL2TZ(f)',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/linnyr419_link_2/LANL2TZF.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']),  # noqa
    'aug-cc-pwCVTZ-DK': ESMLBasisSet(
        'aug-cc-pwCVTZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCVTZ-DK-AGG.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'cc-pCVQZ-F12 MP2 Fitting': ESMLBasisSet(
        'cc-pCVQZ-F12 MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_54/CC-PCVQZ-F12MP2FITTING.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pcJ-1': ESMLBasisSet(
        'aug-pcJ-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_11/AUG-PCJ-1.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'Def2-SVPD': ESMLBasisSet(
        'Def2-SVPD',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/drappoport_new/DEF2-SVPD.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'SARC2-QZV-DKH': ESMLBasisSet(
        'SARC2-QZV-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new_0/SARC2-QZV-DKH.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
    'cc-pV(T+d)Z': ESMLBasisSet(
        'cc-pV(T+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_PLUS_D.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV5Z-NR MP2 Fitting': ESMLBasisSet(
        'cc-pV5Z-NR MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_46/CC-PV5Z-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    '4-22GSP': ESMLBasisSet(
        '4-22GSP',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/4-22GSP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pVTZ MP2 Fitting': ESMLBasisSet(
        'aug-cc-pVTZ MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new/AUG-CC-PVTZ-NRMP2FITTING.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pwCVTZ-DK3': ESMLBasisSet(
        'cc-pwCVTZ-DK3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_37/CC-PWCVTZ-DK3.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']),  # noqa
    'cc-pwCVQZ-DK': ESMLBasisSet(
        'cc-pwCVQZ-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCVQZ-DK.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'maug-cc-pVTZ': ESMLBasisSet(
        'maug-cc-pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Papajak_new_5/MAUG-CC-PVTZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV5Z(fi/sf/lc)': ESMLBasisSet(
        'cc-pV5Z(fi/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z_FI_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pV6Z': ESMLBasisSet(
        'aug-cc-pV6Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PV6Z-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pCV5Z': ESMLBasisSet(
        'aug-cc-pCV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PCV5Z-AGG.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pVDZ': ESMLBasisSet(
        'aug-cc-pVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVDZ-AGG.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVDZ(fi/sf/fw)': ESMLBasisSet(
        'cc-pVDZ(fi/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_FI_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVDZ-F12 MP2 Fitting': ESMLBasisSet(
        'cc-pVDZ-F12 MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_49/CC-PVDZ-F12MP2FITTING.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV5Z-DK': ESMLBasisSet(
        'cc-pV5Z-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z-DK.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Partridge Uncontracted 4': ESMLBasisSet(
        'Partridge Uncontracted 4',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/PARTRIDGE4.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'ccemd-2': ESMLBasisSet(
        'ccemd-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_6/CCEMD-2.xml',
        ['He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar', 'P', 'Cl', 'Al', 'Si', 'S']),  # noqa
    'aug-pcS-0': ESMLBasisSet(
        'aug-pcS-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_19/AUG-PCS-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pV(T+d)Z': ESMLBasisSet(
        'aug-cc-pV(T+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVTPDZ-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pc-2': ESMLBasisSet(
        'aug-pc-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_6/AUG-PC-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'STO-3G*': ESMLBasisSet(
        'STO-3G*',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/STO-3GS-AGG.xml',
        ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCVQZ-RI': ESMLBasisSet(
        'cc-pwCVQZ-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_30/CC-PWCVQZ-RI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    'cc-pVDZ-PP MP2 Fitting': ESMLBasisSet(
        'cc-pVDZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_2/CC-PVDZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'cc-pwCV5Z-PP MP2 Fitting': ESMLBasisSet(
        'cc-pwCV5Z-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_25/CC-PWCV5Z-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'cc-pwCV5Z-RI': ESMLBasisSet(
        'cc-pwCV5Z-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_31/CC-PWCV5Z-RI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    '6-311+G(2d,p)': ESMLBasisSet(
        '6-311+G(2d,p)',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_new_5/6-311PG2DP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'cc-pVQZ': ESMLBasisSet(
        'cc-pVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pCVQZ': ESMLBasisSet(
        'aug-cc-pCVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PCVQZ-AGG.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'SARC-ZORA': ESMLBasisSet(
        'SARC-ZORA',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/demchem_addElements_0/SARC-ZORA.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']),  # noqa
    'aug-cc-pwCVDZ': ESMLBasisSet(
        'aug-cc-pwCVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCVDZ-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVTZ(pt/sf/lc)': ESMLBasisSet(
        'cc-pVTZ(pt/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_PT_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'NASA Ames cc-pCVTZ': ESMLBasisSet(
        'NASA Ames cc-pCVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_CC-PCVTZ.xml',
        ['Ti']),  # noqa
    'pV7Z': ESMLBasisSet(
        'pV7Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/PV7Z.xml',
        ['H', 'C', 'N', 'O', 'F', 'Ne', 'S']),  # noqa
    'aug-pc-0': ESMLBasisSet(
        'aug-pc-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_8/AUG-PC-0.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'QZP-DKH': ESMLBasisSet(
        'QZP-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_21/QZP-DKH.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']),  # noqa
    'pcJ-1_2006': ESMLBasisSet(
        'pcJ-1_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_10/PCJ-1.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV5Z(fi/sf/fw)': ESMLBasisSet(
        'cc-pV5Z(fi/sf/fw)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z_FI_SF_FW.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'ADZP': ESMLBasisSet(
        'ADZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_18/ADZP.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-pcSseg-2': ESMLBasisSet(
        'aug-pcSseg-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_46/AUG-PCSSEG-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '2ZaPa-NR_CV': ESMLBasisSet(
        '2ZaPa-NR_CV',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_21/2ZAPA-NR_CV.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcseg-1': ESMLBasisSet(
        'pcseg-1',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_30/PCSEG-1.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'may-cc-pV(Q+d)Z': ESMLBasisSet(
        'may-cc-pV(Q+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_14/MAY-CC-PVQPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV8Z': ESMLBasisSet(
        'cc-pV8Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV8Z.xml',
        ['H', 'Ne']),  # noqa
    'jun-cc-pV(D+d)Z': ESMLBasisSet(
        'jun-cc-pV(D+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_21/JUN-CC-PVDPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'jul-cc-pV(Q+d)Z': ESMLBasisSet(
        'jul-cc-pV(Q+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_11/JUL-CC-PVQPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-mcc-pVTZ': ESMLBasisSet(
        'aug-mcc-pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-MCC-PVTZ.xml',
        ['H']),  # noqa
    'TZ (Dunning)': ESMLBasisSet(
        'TZ (Dunning)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DUNNINGTZ.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'aug-cc-pwCV5Z-DK': ESMLBasisSet(
        'aug-cc-pwCV5Z-DK',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCV5Z-DK-AGG.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'aug-pcseg-2': ESMLBasisSet(
        'aug-pcseg-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_36/AUG-PCSEG-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Ahlrichs VTZ': ESMLBasisSet(
        'Ahlrichs VTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AHLRICHS_VTZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-pcJ-4': ESMLBasisSet(
        'aug-pcJ-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_14/AUG-PCJ-4.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pwCV5Z-RI tight': ESMLBasisSet(
        'cc-pwCV5Z-RI tight',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_18/CC-PWCV5Z-RITIGHT.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-pcJ-1_2006': ESMLBasisSet(
        'aug-pcJ-1_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_15/AUG-PCJ-1.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pCVTZ-F12': ESMLBasisSet(
        'cc-pCVTZ-F12',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_12/CC-PCVTZ-F12.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-cc-pwCVQZ-PP_OPTRI': ESMLBasisSet(
        'aug-cc-pwCVQZ-PP_OPTRI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_32/AUG-CC-PWCVQZ-PP_OPTRI.xml',
        ['Cu', 'Zn', 'Ag', 'Cd', 'Au', 'Hg']),  # noqa
    '6-311++G': ESMLBasisSet(
        '6-311++G',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/rrenslow_new_3/6-311PPG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'aug-pc-4': ESMLBasisSet(
        'aug-pc-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_4/AUG-PC-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'aug-cc-pwCVQZ-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pwCVQZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_41/AUG-CC-PWCVQZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd']),  # noqa
    'pcS-3': ESMLBasisSet(
        'pcS-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_27/PCS-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVDZ(fi/sf/lc)': ESMLBasisSet(
        'cc-pVDZ(fi/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_FI_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'LANL08(f)': ESMLBasisSet(
        'LANL08(f)',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/linnyr419_link_1/LANL08F.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']),  # noqa
    'TZP': ESMLBasisSet(
        'TZP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_16/TZP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']),  # noqa
    'pcJ-0': ESMLBasisSet(
        'pcJ-0',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_5/PCJ-0.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'STO-6G': ESMLBasisSet(
        'STO-6G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/STO-6G.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'ccJ-pVDZ': ESMLBasisSet(
        'ccJ-pVDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/beud_new/CCJ-PVDZ.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne']),  # noqa
    'MINI (Huzinaga)': ESMLBasisSet(
        'MINI (Huzinaga)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/MINI.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'DZ (Dunning)': ESMLBasisSet(
        'DZ (Dunning)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DUNNINGDZ.xml',
        ['H', 'Li', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    'STO-3G': ESMLBasisSet(
        'STO-3G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/STO-3G.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I']),  # noqa
    'Binning/Curtiss SV': ESMLBasisSet(
        'Binning/Curtiss SV',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/BINCURT-VDZ.xml',
        ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '6-311G': ESMLBasisSet(
        '6-311G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-311G.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'I']),  # noqa
    'cc-pCVDZ-F12': ESMLBasisSet(
        'cc-pCVDZ-F12',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/kipeters_new_11/CC-PCVDZ-F12.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVQZ(fi/sf/sc)': ESMLBasisSet(
        'cc-pVQZ(fi/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVQZ_FI_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pwCV5Z-NR': ESMLBasisSet(
        'cc-pwCV5Z-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PWCV5Z-NR.xml',
        ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']),  # noqa
    'cc-pV5Z': ESMLBasisSet(
        'cc-pV5Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'G3MP2large': ESMLBasisSet(
        'G3MP2large',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/drhaney_new/G3MP2LARGE.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'x2c-SV(P)all-2c': ESMLBasisSet(
        'x2c-SV(P)all-2c',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_3/X2C-SVPALL-2C.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'DZ + Double Rydberg (Dunning-Hay)': ESMLBasisSet(
        'DZ + Double Rydberg (Dunning-Hay)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/DZPDOUBLERYDBERGDUNNING-HAY-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl']),  # noqa
    'pcSseg-3': ESMLBasisSet(
        'pcSseg-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_42/PCSSEG-3.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'NASA Ames cc-pVTZ': ESMLBasisSet(
        'NASA Ames cc-pVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_CC-PVTZ.xml',
        ['Ti', 'Fe']),  # noqa
    'aug-cc-pVQZ-PP MP2 Fitting': ESMLBasisSet(
        'aug-cc-pVQZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_8/AUG-CC-PVQZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'Partridge Uncontracted 1': ESMLBasisSet(
        'Partridge Uncontracted 1',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/PARTRIDGE1.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr']),  # noqa
    'cc-pVDZ(fi/sf/sc)': ESMLBasisSet(
        'cc-pVDZ(fi/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVDZ_FI_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'G3LargeXP': ESMLBasisSet(
        'G3LargeXP',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_24/G3LARGEXP.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'jun-cc-pV(Q+d)Z': ESMLBasisSet(
        'jun-cc-pV(Q+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_13/JUN-CC-PVQPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'pcSseg-2': ESMLBasisSet(
        'pcSseg-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_41/PCSSEG-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'SARC2-QZV-ZORA': ESMLBasisSet(
        'SARC2-QZV-ZORA',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new_2/SARC2-QZV-ZORA.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
    'aug-cc-pV(D+d)Z': ESMLBasisSet(
        'aug-cc-pV(D+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PVDPDZ-AGG.xml',
        ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVTZ-RI': ESMLBasisSet(
        'cc-pVTZ-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_3/CC-PVTZ-RI.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Ahlrichs VDZ': ESMLBasisSet(
        'Ahlrichs VDZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AHLRICHS_VDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'Partridge Uncontracted 2': ESMLBasisSet(
        'Partridge Uncontracted 2',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/PARTRIDGE2.xml',
        ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVTZ(fi/sf/sc)': ESMLBasisSet(
        'cc-pVTZ(fi/sf/sc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PVTZ_FI_SF_SC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVDZ-PP-RI': ESMLBasisSet(
        'cc-pVDZ-PP-RI',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_21/CC-PVDZ-PP-RI.xml',
        ['Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'aug-cc-pVQZ-RI diffuse': ESMLBasisSet(
        'aug-cc-pVQZ-RI diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_9/AUG-CC-PVQZ-RIDIFFUSE.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'SARC2-QZVP-ZORA/JK': ESMLBasisSet(
        'SARC2-QZVP-ZORA/JK',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new_5/SARC2-QZVP-ZORAJK.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
    'apr-cc-pV(Q+d)Z': ESMLBasisSet(
        'apr-cc-pV(Q+d)Z',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/papajak_new_15/APR-CC-PVQPDZ.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pVTZ-PP MP2 Fitting': ESMLBasisSet(
        'cc-pVTZ-PP MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_3/CC-PVTZ-PPMP2FITTING.xml',
        ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt']),  # noqa
    'x2c-SVPall-2c': ESMLBasisSet(
        'x2c-SVPall-2c',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/patrikpollak2704@gmail.com_new_4/X2C-SVPALL-2C.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']),  # noqa
    'Lanl2DZ+2s2p2d2f': ESMLBasisSet(
        'Lanl2DZ+2s2p2d2f',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/omr_link_0/LANL2DZP2S2P2D2F.xml',
        ['Rh']),  # noqa
    'pcJ-3_2006': ESMLBasisSet(
        'pcJ-3_2006',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_12/PCJ-3.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    '6-311++G**-J': ESMLBasisSet(
        '6-311++G**-J',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/sauer_new_3/6-311PPGSS-J.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    'aug-cc-pVTZ-RI diffuse': ESMLBasisSet(
        'aug-cc-pVTZ-RI diffuse',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/haettig_new_8/AUG-CC-PVTZ-RIDIFFUSE.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    '6-31+G*-J': ESMLBasisSet(
        '6-31+G*-J',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/sauer_new_4/6-31PGS-J.xml',
        ['H', 'C', 'N', 'O']),  # noqa
    'cc-pV(T+d)Z+': ESMLBasisSet(
        'cc-pV(T+d)Z+',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Papajak_new_7/CC-PVTPDZP.xml',
        ['H', 'C', 'N', 'O', 'S', 'Cl']),  # noqa
    'aug-mcc-pV6Z': ESMLBasisSet(
        'aug-mcc-pV6Z',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-MCC-PV6Z.xml',
        ['H']),  # noqa
    'SARC2-QZVP-DKH': ESMLBasisSet(
        'SARC2-QZVP-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new/SARC2-QZVP-DKH.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
    'pc-2': ESMLBasisSet(
        'pc-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_correction_2/PC-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'NASA Ames ANO2': ESMLBasisSet(
        'NASA Ames ANO2',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/NASA_ANO2.xml',
        ['Sc', 'Ti', 'V']),  # noqa
    'pSBKJC': ESMLBasisSet(
        'pSBKJC',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/lnvidal_link/PSBKJC.xml',
        ['C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Ge', 'As', 'Se', 'Br', 'Sn', 'Sb', 'Te', 'I']),  # noqa
    '4ZaPa-NR': ESMLBasisSet(
        '4ZaPa-NR',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/comartin_new_4/4ZAPA-NR.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'aug-mcc-pVQZ': ESMLBasisSet(
        'aug-mcc-pVQZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-MCC-PVQZ.xml',
        ['H']),  # noqa
    'aug-pcSseg-4': ESMLBasisSet(
        'aug-pcSseg-4',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/frj_new_48/AUG-PCSSEG-4.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'cc-pVQZ-F12 MP2 Fitting': ESMLBasisSet(
        'cc-pVQZ-F12 MP2 Fitting',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/GHill_new_51/CC-PVQZ-F12MP2FITTING.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'coemd-2': ESMLBasisSet(
        'coemd-2',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_2/COEMD-2.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'MIDI!': ESMLBasisSet(
        'MIDI!',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/MIDI_BANG.xml',
        ['H', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I']),  # noqa
    '6ZP-DKH': ESMLBasisSet(
        '6ZP-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_23/6ZP-DKH.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'ccemd-3': ESMLBasisSet(
        'ccemd-3',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/jussilehtola_new_7/CCEMD-3.xml',
        ['He', 'H', 'Be', 'Li', 'Ne', 'C', 'N', 'F', 'B', 'O', 'Mg', 'Na', 'Ar', 'P', 'Cl', 'Al', 'Si', 'S']),  # noqa
    '6-311++G**': ESMLBasisSet(
        '6-311++G**',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/6-311PPGSS-AGG.xml',
        ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']),  # noqa
    'aug-cc-pwCVTZ': ESMLBasisSet(
        'aug-cc-pwCVTZ',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/AUG-CC-PWCVTZ-AGG.xml',
        ['B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']),  # noqa
    'cc-pV5Z(pt/sf/lc)': ESMLBasisSet(
        'cc-pV5Z(pt/sf/lc)',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/CC-PV5Z_PT_SF_LC.xml',
        ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']),  # noqa
    'SARC-DKH': ESMLBasisSet(
        'SARC-DKH',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/ganymede_new/SARC-DKH.xml',
        ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']),  # noqa
    '4-31G': ESMLBasisSet(
        '4-31G',
        '/files/projects/Basis_Set_Curators/Gaussian/emsl-lib/4-31G.xml',
        ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'P', 'S', 'Cl']),  # noqa
    'SARC2-QZVP-ZORA': ESMLBasisSet(
        'SARC2-QZVP-ZORA',
        '/files/projects/Basis_Set_Curators/Gaussian/contrib/Daniel Aravena_new_1/SARC2-QZVP-ZORA.xml',
        ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']),  # noqa
}

# and and ordered list of the keys
AVAILABLE_BS_NAME = [
    '2ZaPa-NR',
    '2ZaPa-NR_CV',
    '3-21++G',
    '3-21++G*',
    '3-21G',
    '3-21G*',
    '3-21GSP',
    '3ZaPa-NR',
    '3ZaPa-NR_CV',
    '4-22GSP',
    '4-31G',
    '4ZaPa-NR',
    '4ZaPa-NR_CV',
    '5ZP',
    '5ZP-DKH',
    '5ZaPa-NR',
    '5ZaPa-NR_CV',
    '6-31++G',
    '6-31++G*',
    '6-31++G**',
    '6-31++G**-J',
    '6-31+G',
    '6-31+G*',
    '6-31+G**',
    '6-31+G*-J',
    '6-311++G',
    '6-311++G(2d,2p)',
    '6-311++G(3df,3pd)',
    '6-311++G*',
    '6-311++G**',
    '6-311++G**-J',
    '6-311+G',
    '6-311+G(2d,p)',
    '6-311+G*',
    '6-311+G**',
    '6-311+G*-J',
    '6-311G',
    '6-311G(2df,2pd)',
    '6-311G*',
    '6-311G**',
    '6-311G-J',
    '6-31G',
    '6-31G(2df,p)',
    '6-31G(3df,3pd)',
    '6-31G*',
    '6-31G**',
    '6-31G*-Blaudeau',
    '6-31G-Blaudeau',
    '6-31G-J',
    '6ZP',
    '6ZP-DKH',
    '6ZaPa-NR',
    '7ZaPa-NR',
    'A5ZP',
    'ADZP',
    'ANO-RCC',
    'AQZP',
    'ATZP',
    'Ahlrichs TZV',
    'Ahlrichs VDZ',
    'Ahlrichs VTZ',
    'Ahlrichs pVDZ',
    'B2 basis set for Zn',
    'Bauschlicher ANO',
    'Binning/Curtiss SV',
    'Binning/Curtiss SVP',
    'Binning/Curtiss VTZ',
    'Binning/Curtiss VTZP',
    'CVTZ',
    'Chipman DZP + Diffuse',
    'Cologne DKH2',
    'DZ (Dunning)',
    'DZ + Double Rydberg (Dunning-Hay)',
    'DZ + Rydberg (Dunning)',
    'DZP',
    'DZP (Dunning)',
    'DZP + Diffuse (Dunning)',
    'DZP + Rydberg (Dunning)',
    'DZP-DKH',
    'DZQ',
    'Def2-QZVP',
    'Def2-QZVPD',
    'Def2-QZVPP',
    'Def2-QZVPPD',
    'Def2-SV(P)',
    'Def2-SVP',
    'Def2-SVPD',
    'Def2-TZVP',
    'Def2-TZVPD',
    'Def2-TZVPP',
    'Def2-TZVPPD',
    'Feller Misc. CVDZ',
    'Feller Misc. CVQZ',
    'Feller Misc. CVTZ',
    'G3LargeXP',
    'G3MP2LargeXP',
    'G3MP2large',
    'GAMESS PVTZ',
    'GAMESS VTZ',
    'IGLO-II',
    'IGLO-III',
    'LANL08',
    'LANL08(f)',
    'LANL08+',
    'LANL08d',
    'LANL2TZ',
    'LANL2TZ(f)',
    'LANL2TZ+',
    'Lanl2-[10s8p7d3f2g]',
    'Lanl2-[5s4p4d2f]',
    'Lanl2-[6s4p4d2f]',
    'Lanl2DZ+1d1f',
    'Lanl2DZ+2s2p2d2f',
    'MG3S',
    'MIDI (Huzinaga)',
    'MIDI!',
    'MINI (Huzinaga)',
    'MINI (Scaled)',
    'MINIs-BSIP1',
    'McLean/Chandler VTZ',
    'NASA Ames ANO',
    'NASA Ames ANO2',
    'NASA Ames cc-pCV5Z',
    'NASA Ames cc-pCVQZ',
    'NASA Ames cc-pCVTZ',
    'NASA Ames cc-pV5Z',
    'NASA Ames cc-pVQZ',
    'NASA Ames cc-pVTZ',
    'NLO-V',
    'NMR-DKH (TZ2P)',
    'Partridge Uncontracted 1',
    'Partridge Uncontracted 2',
    'Partridge Uncontracted 3',
    'Partridge Uncontracted 4',
    'Pt - mDZP',
    'QZP',
    'QZP-DKH',
    'Roos Augmented Double Zeta ANO',
    'Roos Augmented Triple Zeta ANO',
    'SARC-DKH',
    'SARC-ZORA',
    'SARC2-QZV-DKH',
    'SARC2-QZV-DKH/JK',
    'SARC2-QZV-ZORA',
    'SARC2-QZV-ZORA/JK',
    'SARC2-QZVP-DKH',
    'SARC2-QZVP-DKH/JK',
    'SARC2-QZVP-ZORA',
    'SARC2-QZVP-ZORA/JK',
    'SBKJC Polarized (p,2d) - LFK',
    'STO-2G',
    'STO-3G',
    'STO-3G*',
    'STO-6G',
    'SV (Dunning-Hay)',
    'SV + Double Rydberg (Dunning-Hay)',
    'SV + Rydberg (Dunning-Hay)',
    'SVP (Dunning-Hay)',
    'SVP + Diffuse (Dunning-Hay)',
    'SVP + Diffuse + Rydberg',
    'SVP + Rydberg (Dunning-Hay)',
    'Sadlej pVTZ',
    'Sadlej+',
    'TZ (Dunning)',
    'TZP',
    'TZP-DKH',
    'UGBS',
    'WTBS',
    'Wachters+f',
    'Weigend Coulomb Fitting',
    'apr-cc-pV(Q+d)Z',
    'aug-cc-pCV(T+d)Z',
    'aug-cc-pCV5Z',
    'aug-cc-pCVDZ',
    'aug-cc-pCVQZ',
    'aug-cc-pCVTZ',
    'aug-cc-pV(5+d)Z',
    'aug-cc-pV(6+d)Z',
    'aug-cc-pV(7+d)Z',
    'aug-cc-pV(D+d)Z',
    'aug-cc-pV(Q+d)Z',
    'aug-cc-pV(T+d)Z',
    'aug-cc-pV5Z',
    'aug-cc-pV5Z-DK',
    'aug-cc-pV5Z-NR MP2 Fitting',
    'aug-cc-pV5Z-PP MP2 Fitting',
    'aug-cc-pV5Z-PP-RI Diffuse',
    'aug-cc-pV5Z-PP_OPTRI',
    'aug-cc-pV5Z-RI diffuse',
    'aug-cc-pV5Z_OPTRI',
    'aug-cc-pV6Z',
    'aug-cc-pV6Z-RI diffuse',
    'aug-cc-pV7Z',
    'aug-cc-pVDZ',
    'aug-cc-pVDZ-DK',
    'aug-cc-pVDZ-PP MP2 Fitting',
    'aug-cc-pVDZ-PP-RI Diffuse',
    'aug-cc-pVDZ-PP_OPTRI',
    'aug-cc-pVDZ-RI diffuse',
    'aug-cc-pVDZ_OPTRI',
    'aug-cc-pVQZ',
    'aug-cc-pVQZ-DK',
    'aug-cc-pVQZ-NR MP2 Fitting',
    'aug-cc-pVQZ-PP MP2 Fitting',
    'aug-cc-pVQZ-PP-RI Diffuse',
    'aug-cc-pVQZ-PP_OPTRI',
    'aug-cc-pVQZ-RI diffuse',
    'aug-cc-pVQZ_OPTRI',
    'aug-cc-pVTZ',
    'aug-cc-pVTZ MP2 Fitting',
    'aug-cc-pVTZ-DK',
    'aug-cc-pVTZ-J',
    'aug-cc-pVTZ-PP MP2 Fitting',
    'aug-cc-pVTZ-PP-RI Diffuse',
    'aug-cc-pVTZ-PP_OPTRI',
    'aug-cc-pVTZ-RI diffuse',
    'aug-cc-pVTZ_OPTRI',
    'aug-cc-pwCV5Z',
    'aug-cc-pwCV5Z-DK',
    'aug-cc-pwCV5Z-NR',
    'aug-cc-pwCV5Z-PP MP2 Fitting',
    'aug-cc-pwCV5Z-PP_OPTRI',
    'aug-cc-pwCVDZ',
    'aug-cc-pwCVDZ-PP MP2 Fitting',
    'aug-cc-pwCVDZ-PP_OPTRI',
    'aug-cc-pwCVQZ',
    'aug-cc-pwCVQZ-DK',
    'aug-cc-pwCVQZ-NR',
    'aug-cc-pwCVQZ-PP MP2 Fitting',
    'aug-cc-pwCVQZ-PP_OPTRI',
    'aug-cc-pwCVTZ',
    'aug-cc-pwCVTZ-DK',
    'aug-cc-pwCVTZ-NR',
    'aug-cc-pwCVTZ-NR MP2 Fitting',
    'aug-cc-pwCVTZ-PP MP2 Fitting',
    'aug-cc-pwCVTZ-PP_OPTRI',
    'aug-mcc-pV5Z',
    'aug-mcc-pV6Z',
    'aug-mcc-pV7Z',
    'aug-mcc-pV8Z',
    'aug-mcc-pVQZ',
    'aug-mcc-pVTZ',
    'aug-pV7Z',
    'aug-pc-0',
    'aug-pc-1',
    'aug-pc-2',
    'aug-pc-3',
    'aug-pc-4',
    'aug-pcJ-0',
    'aug-pcJ-0_2006',
    'aug-pcJ-1',
    'aug-pcJ-1_2006',
    'aug-pcJ-2',
    'aug-pcJ-2_2006',
    'aug-pcJ-3',
    'aug-pcJ-3_2006',
    'aug-pcJ-4',
    'aug-pcJ-4_2006',
    'aug-pcS-0',
    'aug-pcS-1',
    'aug-pcS-2',
    'aug-pcS-3',
    'aug-pcS-4',
    'aug-pcSseg-0',
    'aug-pcSseg-1',
    'aug-pcSseg-2',
    'aug-pcSseg-3',
    'aug-pcSseg-4',
    'aug-pcseg-0',
    'aug-pcseg-1',
    'aug-pcseg-2',
    'aug-pcseg-3',
    'aug-pcseg-4',
    'cc-pCV5Z',
    'cc-pCV6Z',
    'cc-pCV6Z (old)',
    'cc-pCV6Z(old)',
    'cc-pCVDZ',
    'cc-pCVDZ(old)',
    'cc-pCVDZ-F12',
    'cc-pCVDZ-F12 MP2 Fitting',
    'cc-pCVDZ-F12_OPTRI',
    'cc-pCVQZ',
    'cc-pCVQZ(old)',
    'cc-pCVQZ-F12',
    'cc-pCVQZ-F12 MP2 Fitting',
    'cc-pCVQZ-F12_OPTRI',
    'cc-pCVTZ',
    'cc-pCVTZ(old)',
    'cc-pCVTZ-F12',
    'cc-pCVTZ-F12 MP2 Fitting',
    'cc-pCVTZ-F12_OPTRI',
    'cc-pV(5+d)Z',
    'cc-pV(5+d)Z-RI',
    'cc-pV(6+d)Z',
    'cc-pV(6+d)Z-RI',
    'cc-pV(D+d)Z',
    'cc-pV(D+d)Z-RI',
    'cc-pV(Q+d)Z',
    'cc-pV(Q+d)Z-RI',
    'cc-pV(T+d)Z',
    'cc-pV(T+d)Z+',
    'cc-pV(T+d)Z-DK',
    'cc-pV(T+d)Z-RI',
    'cc-pV5Z',
    'cc-pV5Z(fi/sf/fw)',
    'cc-pV5Z(fi/sf/lc)',
    'cc-pV5Z(fi/sf/sc)',
    'cc-pV5Z(pt/sf/fw)',
    'cc-pV5Z(pt/sf/lc)',
    'cc-pV5Z(pt/sf/sc)',
    'cc-pV5Z-DK',
    'cc-pV5Z-NR MP2 Fitting',
    'cc-pV5Z-PP MP2 Fitting',
    'cc-pV5Z-PP-RI',
    'cc-pV5Z-RI',
    'cc-pV6Z',
    'cc-pV6Z-RI',
    'cc-pV8Z',
    'cc-pV9Z',
    'cc-pVDZ',
    'cc-pVDZ(fi/sf/fw)',
    'cc-pVDZ(fi/sf/lc)',
    'cc-pVDZ(fi/sf/sc)',
    'cc-pVDZ(pt/sf/fw)',
    'cc-pVDZ(pt/sf/lc)',
    'cc-pVDZ(pt/sf/sc)',
    'cc-pVDZ(seg-opt)',
    'cc-pVDZ-DK',
    'cc-pVDZ-DK3',
    'cc-pVDZ-F12',
    'cc-pVDZ-F12 MP2 Fitting',
    'cc-pVDZ-F12_OPTRI',
    'cc-pVDZ-PP MP2 Fitting',
    'cc-pVDZ-PP-RI',
    'cc-pVDZ-RI',
    'cc-pVQZ',
    'cc-pVQZ(fi/sf/fw)',
    'cc-pVQZ(fi/sf/lc)',
    'cc-pVQZ(fi/sf/sc)',
    'cc-pVQZ(pt/sf/fw)',
    'cc-pVQZ(pt/sf/lc)',
    'cc-pVQZ(pt/sf/sc)',
    'cc-pVQZ(seg-opt)',
    'cc-pVQZ-DK',
    'cc-pVQZ-DK3',
    'cc-pVQZ-F12',
    'cc-pVQZ-F12 MP2 Fitting',
    'cc-pVQZ-F12_OPTRI',
    'cc-pVQZ-NR MP2 Fitting',
    'cc-pVQZ-PP MP2 Fitting',
    'cc-pVQZ-PP-RI',
    'cc-pVQZ-RI',
    'cc-pVTZ',
    'cc-pVTZ MP2 Fitting',
    'cc-pVTZ(fi/sf/fw)',
    'cc-pVTZ(fi/sf/lc)',
    'cc-pVTZ(fi/sf/sc)',
    'cc-pVTZ(pt/sf/fw)',
    'cc-pVTZ(pt/sf/lc)',
    'cc-pVTZ(pt/sf/sc)',
    'cc-pVTZ(seg-opt)',
    'cc-pVTZ-DK',
    'cc-pVTZ-DK3',
    'cc-pVTZ-F12',
    'cc-pVTZ-F12 MP2 Fitting',
    'cc-pVTZ-F12_OPTRI',
    'cc-pVTZ-PP MP2 Fitting',
    'cc-pVTZ-PP-RI',
    'cc-pVTZ-RI',
    'cc-pwCV5Z',
    'cc-pwCV5Z Core Set',
    'cc-pwCV5Z-DK',
    'cc-pwCV5Z-NR',
    'cc-pwCV5Z-PP MP2 Fitting',
    'cc-pwCV5Z-RI',
    'cc-pwCV5Z-RI tight',
    'cc-pwCVDZ',
    'cc-pwCVDZ-DK3',
    'cc-pwCVDZ-PP MP2 Fitting',
    'cc-pwCVDZ-RI',
    'cc-pwCVDZ-RI tight',
    'cc-pwCVQZ',
    'cc-pwCVQZ-DK',
    'cc-pwCVQZ-DK3',
    'cc-pwCVQZ-NR',
    'cc-pwCVQZ-PP MP2 Fitting',
    'cc-pwCVQZ-RI',
    'cc-pwCVQZ-RI tight',
    'cc-pwCVTZ',
    'cc-pwCVTZ-DK',
    'cc-pwCVTZ-DK3',
    'cc-pwCVTZ-NR',
    'cc-pwCVTZ-NR MP2 Fitting',
    'cc-pwCVTZ-PP MP2 Fitting',
    'cc-pwCVTZ-RI',
    'cc-pwCVTZ-RI tight',
    'ccJ-pV5Z',
    'ccJ-pVDZ',
    'ccJ-pVQZ',
    'ccJ-pVTZ',
    'ccemd-2',
    'ccemd-3',
    'coemd-2',
    'coemd-3',
    'coemd-4',
    'coemd-ref',
    'd-aug-cc-pV5Z',
    'd-aug-cc-pV6Z',
    'd-aug-cc-pVDZ',
    'd-aug-cc-pVQZ',
    'd-aug-cc-pVTZ',
    'dhf-QZVP',
    'dhf-QZVPP',
    'dhf-SV(P)',
    'dhf-SVP',
    'dhf-TZVP',
    'dhf-TZVPP',
    'jul-cc-pV(D+d)Z',
    'jul-cc-pV(Q+d)Z',
    'jul-cc-pV(T+d)Z',
    'jun-cc-pV(D+d)Z',
    'jun-cc-pV(Q+d)Z',
    'jun-cc-pV(T+d)Z',
    'm6-31G',
    'm6-31G*',
    'maug-cc-pV(D+d)Z',
    'maug-cc-pV(Q+d)Z',
    'maug-cc-pV(T+d)Z',
    'maug-cc-pVDZ',
    'maug-cc-pVQZ',
    'maug-cc-pVTZ',
    'may-cc-pV(Q+d)Z',
    'may-cc-pV(T+d)Z',
    'modified LANL2DZ',
    'pSBKJC',
    'pV6Z',
    'pV7Z',
    'pc-0',
    'pc-1',
    'pc-2',
    'pc-3',
    'pc-4',
    'pcJ-0',
    'pcJ-0_2006',
    'pcJ-1',
    'pcJ-1_2006',
    'pcJ-2',
    'pcJ-2_2006',
    'pcJ-3',
    'pcJ-3_2006',
    'pcJ-4',
    'pcJ-4_2006',
    'pcS-0',
    'pcS-1',
    'pcS-2',
    'pcS-3',
    'pcS-4',
    'pcSseg-0',
    'pcSseg-1',
    'pcSseg-2',
    'pcSseg-3',
    'pcSseg-4',
    'pcemd-2',
    'pcemd-3',
    'pcemd-4',
    'pcseg-0',
    'pcseg-1',
    'pcseg-2',
    'pcseg-3',
    'pcseg-4',
    's3-21G',
    's3-21G*',
    's6-31G',
    's6-31G*',
    'un-ccemd-ref',
    'un-pcemd-ref',
    'x2c-SV(P)all',
    'x2c-SV(P)all-2c',
    'x2c-SVPall',
    'x2c-SVPall-2c',
    'x2c-TZVPPall',
    'x2c-TZVPPall-2c',
    'x2c-TZVPall',
    'x2c-TZVPall-2c',
    'x2c-coulomb-fitting'
]
