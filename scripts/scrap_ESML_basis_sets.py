#!/usr/bin/env python3
"""
Update the list of basis sets from ESML basis set exchange.
"""

import argparse
import requests
import os
from bs4 import BeautifulSoup
from datetime import datetime

JS_PEID = 11535052407933
SCRAP_URL = 'https://bse.pnl.gov/bse/portal/user/anon/js_peid/{}/panel/Main/template/content'
DESTINATION_FILE = 'qcip_tools/basis_set_esml.py'
DEFINITION_BS = 'new basisSet('

FILE_CONTENT = """{file_beginning}# DO NOT EDIT BELOW THIS POINT!
# generated on {date}


#: current JS PEID
JS_PEID = {jspeid}

#: URL for the download
URL = 'https://bse.pnl.gov:443/bse/portal/user/anon/js_peid/{{}}/action/' \\
      'portlets.BasisSetAction/template/courier_content/panel/Main/eventSubmit_doDownload/true'

#: List of authorized formats
AUTHORIZED_FORMATS = [
{authorized_formats}]

# the basis sets:
AVAILABLE_BS = {{
{dict}}}

# and and ordered list of the keys
AVAILABLE_BS_NAME = [
{list}]
"""

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


parser = argparse.ArgumentParser(__doc__)

parser.add_argument(
    '-i', '--jspeid', action='store', default=JS_PEID, help='js peid')

parser.add_argument(
    '-d', '--destination', action='store', default=DESTINATION_FILE, help='where to write the basis sets')


def main():
    args = parser.parse_args()
    print(__doc__)

    if not os.path.exists(DESTINATION_FILE):
        print('cannot access {}'.format(args.destination))
        return False

    url = SCRAP_URL.format(args.jspeid)

    print('Requesting (may be long) ... ')

    r = requests.get(url)
    if r.status_code != 200:
        print('Connection status code is {} (try with a different JS_PEID?)'.format(r.status_code))
        return False

    print('... ok!')

    soup = BeautifulSoup(r.content, 'html.parser')
    script_parts = soup.find_all('script', {'language': 'JavaScript'})

    basis_sets = {}

    for script in script_parts:
        s = script.string
        if 'function basisSet(' in s:
            found_bs = s.find(DEFINITION_BS, 0)
            while found_bs != -1:
                found_end = s.find(');', found_bs)
                if found_end == -1:
                    raise Exception('found beginning but no end ?')
                found_lst = s.find('"[', found_bs)
                if found_lst == -1:
                    raise Exception('cannot find list of atom?')
                found_lst_end = s.find(']"', found_lst)
                if found_lst_end == -1:
                    raise Exception('cannot find end of list of atoms')

                definition = s[found_bs + len(DEFINITION_BS):found_end]
                relative_lst_end = found_lst_end-found_bs-len(DEFINITION_BS)
                x = definition[:relative_lst_end].split('"')
                y = definition[relative_lst_end:].split(',')
                path, name, bs_type, atoms = x[1], x[3], x[5], x[7][1:]

                status = y[1][2:-1]

                if bs_type == 'orbital' and status == 'published':
                    if name not in basis_sets:
                        basis_sets[name] = [path, atoms]
                    else:
                        print('!! duplicate {}'.format(name))

                found_bs = s.find(DEFINITION_BS, found_bs + 1)
            break

    if not basis_sets:
        raise Exception('failed to find basis sets in the ESML page (try with a different JS_PEID?)')

    # catch all authorized formats
    authorized_formats = []

    def catch_option(tag):
        return tag.name == 'option' and \
               tag.parent.name == 'select' and \
               tag.parent.has_attr('name') and \
               tag.parent['name'] == 'outputcode'

    for x in soup.find_all(catch_option):
        authorized_formats.append(x.string)

    # catch file beginning
    file_beginning = ''

    try:
        with open(DESTINATION_FILE) as f:
            while True:
                l = f.readline()
                if '# DO NOT EDIT BELOW THIS POINT!' in l:
                    break
                file_beginning += l
    except IOError as e:
        print(e)
        return False

    # write all stuffs
    try:
        with open(DESTINATION_FILE, 'w') as f:

            content_dict = ''
            content_list = []

            for name, (path, atoms) in basis_sets.items():
                content_dict += \
                    '    \'{n}\': ESMLBasisSet(\n' \
                    '        \'{n}\',\n' \
                    '        \'{p}\',\n' \
                    '        {a}\n'.format(
                        n=name,
                        p=path,
                        a='[{}]),  # noqa'.format(', '.join('\'' + x.strip() + '\'' for x in atoms.split(','))))

                content_list.append(name)

            content_list.sort()

            f.write(
                FILE_CONTENT.format(
                    file_beginning=file_beginning,
                    date=datetime.now().strftime('%B %d, %Y (%H:%M)'),
                    jspeid=args.jspeid,
                    dict=content_dict,
                    list='    \'' + ('\',\n    \''.join(content_list)) + '\'\n',
                    authorized_formats='    \'' + ('\',\n    \''.join(authorized_formats)) + '\'\n'
                ))

            print('\nDone: {} published (!) orbital basis sets were retrieved '
                  '(excluding ECP, polarization, diffuse ...)'.format(len(basis_sets)))
    except IOError as e:
        print(e)
        return False


if __name__ == '__main__':
    main()
