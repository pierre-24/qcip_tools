#!/usr/bin/env python3
"""
Make the release.

Note: make sure that you have the rights on the target branch of both {} and {}.
"""

import argparse
import sys
import subprocess
from datetime import datetime

import qcip_tools

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

SOURCE_BRANCH = 'dev'
TARGET_BRANCH = 'master'
GIT_REMOTE = 'git.pierrebeaujean.net:pierre/'
UNAMUR_REMOTE = 'gitlab.unamur.be:pierre.beaujean/'
STABLE_VERSION = 'Stable version: '
RELEASE_TAG = 'release-v{version}'
URL_VERSION = 'https://git.pierrebeaujean.net/pierre/qcip_tools/tree/{}'.format(RELEASE_TAG)
RELEASE_COMMIT = '[{}] Upgrade the README for the release'


class ValidationError(ValueError):
    pass


def input_with_validation(msg, validators=None, default_value=None):
    """Check input

    :param msg: message to ask
    :type msg: str
    :param validators: validator functions
    :type validators: list
    :param default_value: default value
    :type default_value: str
    :rtype: str
    """

    while True:
        try:
            r = input('{} ? {}'.format(msg, ' [{}] '.format(default_value) if default_value else ''))
        except (KeyboardInterrupt, EOFError):
            print('aborting')
            sys.exit(0)

        if r == '' and default_value:
            r = default_value

        if validators:
            try:
                for validator in validators:
                    r = validator(r)
                break
            except ValidationError as e:
                print('ERROR:', e)
                continue
        else:
            break
    return r


def not_empty_validator(r):
    if not r:
        raise ValidationError('empty input')
    return r


def yesorno_validator(r):
    means_yes = ['y', 'Y', 'yes', 'YES']
    means_no = ['N', 'n', 'no', 'NO']

    if r in means_yes:
        return 'Y'
    elif r in means_no:
        return 'N'
    else:
        raise ValidationError('{} is not yes or no'.format(r))


def run_git_command(parameters, forget_exception=False):
    """Just "git xxx"

    :param parameters: list of parameters (Popen style)
    :type parameters: list
    :param forget_exception: sometimes, git use stderr instead of stdout
    :type forget_exception: bool
    :rtype: str
    """
    cmd = ['git']
    cmd.extend(parameters)

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_t, stderr_t = process.communicate()

    if len(stderr_t) != 0:
        if not forget_exception:
            raise Exception('Error: {}'.format(stderr_t.decode()))
        else:
            return stderr_t.decode()

    return stdout_t.decode()


def get_git_config(conf):
    """Just "git config --get xxx"

    :param conf: config string
    :type conf: str
    :rtype: str
    """
    return run_git_command(['config', '--get', conf])

parser = argparse.ArgumentParser('Make the release')

parser.add_argument(
    '-s', '--source', action='store', default=SOURCE_BRANCH, help='Source branch')
parser.add_argument(
    '-t', '--target', action='store', default=TARGET_BRANCH, help='Source branch')
parser.add_argument(
    '-r', '--release', action='store', default=qcip_tools.__version__, help='Release version')


def main():
    args = parser.parse_args()
    print(__doc__.format(GIT_REMOTE, UNAMUR_REMOTE))

    print('Source branch = {}, target branch = {}, release = {}'.format(
        args.source, args.target, args.release))

    if input_with_validation(
            'Are those parameters ok?', validators=[yesorno_validator], default_value='Y') == 'N':
        return False

    tag = RELEASE_TAG.format(version=args.release)

    print('\nOk, lets do the release on "{}":\n'.format(tag))

    # find remote
    remote_git = ''
    remote_unamur = ''

    print('- find remotes ', end='')

    remotes = run_git_command(['remote'])
    for remote in remotes.split():
        remote_url = get_git_config('remote.{}.url'.format(remote))
        if GIT_REMOTE in remote_url:
            remote_git = remote
        if UNAMUR_REMOTE in remote_url:
            remote_unamur = remote
    if not remote_unamur:
        print('[NOK]: cannot find remote for {}'.format(UNAMUR_REMOTE))
        return False
    if not remote_git:
        print('[NOK]: cannot find remote for {}'.format(GIT_REMOTE))
        return False

    print('[OK]: "{}" and "{}" (= UNamur\'s gitlab)'.format(remote_git, remote_unamur))

    # check track on target branch
    print('- check that source ({}) branch tracks {} '.format(args.source, remote_git), end='')
    try:
        target_remote = get_git_config('branch.{}.remote'.format(args.source))
    except Exception as e:
        print('[NOK]: {}'.format(e))
        return False

    if target_remote.strip() != remote_git:
        print('[NOK]: it tracks {}'.format(target_remote))
    else:
        print('[OK]')

    print('- check that target ({}) branch tracks {} '.format(args.target, remote_git), end='')
    try:
        target_remote = get_git_config('branch.{}.remote'.format(args.target))
    except Exception as e:
        print('[NOK]: {}'.format(e))
        return False

    if target_remote.strip() != remote_git:
        print('[NOK]: it tracks {}'.format(target_remote))
    else:
        print('[OK]')

    # fetching
    print('- fetching ', end='')
    run_git_command(['fetch', remote_git], forget_exception=True)
    print('[OK]')

    # checkout to target
    current_branch = run_git_command(['rev-parse', '--abbrev-ref', 'HEAD']).strip()

    if current_branch != args.source:
        print('- checkout source ', end='')
        try:
            run_git_command(['checkout', args.source])
        except Exception as e:
            if 'Switched' not in str(e):
                print('[NOK]: {}'.format(e))
                return False
        print('[OK]')

    # README.md
    print('- pulling on {} '.format(args.source), end='')
    run_git_command(['pull'], forget_exception=True)
    print('[OK]')

    print('- upgrading README.md ', end='')

    readme_content = ''
    try:
        with open('README.md') as f:
            readme_content = f.readlines()
    except IOError as e:
        print('[NOK]: {}'.format(e))

    previous_release = ''

    for index, l in enumerate(readme_content):
        if '<!-- STABLE: -->' in l:
            previous_release = readme_content[index + 1].replace(STABLE_VERSION, '')
            now = datetime.now()
            readme_content[index + 1] = '{}[`{}`]({}) ({})\n'.format(
                STABLE_VERSION,
                tag,
                URL_VERSION.format(version=args.release),
                now.strftime('%B %d, %Y'))
            continue
        if '<!-- PREVIOUS: -->' in l:
            if not previous_release:
                print('[NOK]: previous release (<!-- PREVIOUS: -->) not found in README.md')
                return False
            readme_content.insert(index + 1, '+ {}'.format(previous_release))
            break

    try:
        with open('README.md', 'w') as f:
            f.write(''.join(readme_content))
    except IOError as e:
        print('[NOK]: {}'.format(e))
        return False
    print('[OK]')

    print('- commiting ', end='')
    try:
        run_git_command(['commit', '-am', RELEASE_COMMIT.format(tag)])
    except Exception as e:
        print('[NOK]: {}'.format(e))
        return False
    print('[OK]')

    # checkout target and merge
    print('- checkout target ', end='')
    try:
        run_git_command(['checkout', args.target])
    except Exception as e:
        if 'Switched' not in str(e):
            print('[NOK]: {}'.format(e))
            return False
    print('[OK]')

    print('- merging source ({}) into target ({}) '.format(args.source, args.target), end='')
    try:
        run_git_command(['merge', '--no-ff', args.source])
    except Exception as e:
        print('[NOK]: {}'.format(e))
        return False
    print('[OK]')

    # tagging and pushing!
    print('- tagging ', end='')
    try:
        last_commit = run_git_command(['rev-parse', 'HEAD']).strip()
        run_git_command(['tag', tag, last_commit])
    except Exception as e:
        print('[NOK]: {}'.format(e))
        return False
    print('[OK]: tagging {} on {}'.format(tag, last_commit[:8]))

    print('- pushing ', end='')
    run_git_command(['push', remote_git, args.source], forget_exception=True)
    print('source, ', end='')
    run_git_command(['push', remote_git, args.target], forget_exception=True)
    print('target, ', end='')
    run_git_command(['push', remote_git, tag], forget_exception=True)
    print('tag, ', end='')
    run_git_command(['push', remote_unamur, args.target], forget_exception=True)
    print('target on {}, '.format(remote_unamur), end='')
    run_git_command(['push', remote_unamur, tag], forget_exception=True)
    print('tag on {} '.format(remote_unamur), end='')
    print('[OK]')

    print('\nOK, version {} was released!'.format(args.release))

    # go back to first branch:
    run_git_command(['checkout', current_branch], forget_exception=True)

if __name__ == '__main__':
    main()
