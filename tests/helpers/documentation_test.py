#!/usr/bin/env python3
"""Portable checks for public documentation links and the genotyper manual contract."""
import pathlib
import re
import sys
from urllib.parse import unquote, urlsplit


def prose(text):
    """Ignore fenced examples when extracting actual Markdown links/headings."""
    lines, fenced = [], False
    for line in text.splitlines():
        if line.startswith('```'):
            fenced = not fenced
        elif not fenced:
            lines.append(line)
    assert not fenced, 'Unclosed code fence'
    return '\n'.join(lines)


def headings(text):
    result, seen = set(), {}
    for title in re.findall(r'^#{1,6}\s+(.+?)\s*#*$', prose(text), re.M):
        title = re.sub(r'<[^>]*>', '', title).lower()
        slug = re.sub(r'[^\w\- ]', '', title).replace(' ', '-')
        count = seen.get(slug, 0)
        seen[slug] = count + 1
        result.add(slug + (f'-{count}' if count else ''))
    return result


def check(root):
    pages = [root/'README.md', *sorted((root/'docs').rglob('*.md'))]
    links = 0
    for path in pages:
        text = path.read_text()
        body = prose(text)
        for reference in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', body):
            target = reference.strip().split(' "', 1)[0].strip('<>')
            parsed = urlsplit(target)
            if parsed.scheme or parsed.netloc:
                continue
            destination = (path.parent/unquote(parsed.path)).resolve() if parsed.path else path
            assert destination.exists(), f'{path.relative_to(root)}: missing link target {target}'
            if parsed.fragment and destination.suffix == '.md':
                assert unquote(parsed.fragment) in headings(destination.read_text()), (
                    f'{path.relative_to(root)}: missing heading {target}')
            links += 1
        for retired in ('docs/guides/', 'docs/releases/', '../guides/', '../releases/',
                        'tests/analysis/', 'docs/reports/', '/Users/'):
            assert retired not in text, f'{path.relative_to(root)}: obsolete/private dependency {retired}'
    assert not (root/'docs/guides').exists() and not (root/'docs/releases').exists()
    module = (root/'docs/modules/genotype.md').read_text()
    titles = re.findall(r'^## (.+)$', module, re.M)
    assert titles == ['What it does', 'Required inputs', 'Key options', 'Outputs', 'Limitations', 'Example'], titles
    for term in ('PREDICTION', 'PROVISIONAL', 'NO_PREDICTION', 'GENOTYPE_RESOLVED_SEQUENCE_SET',
                 'SEQUENCE_TIED', 'NOT_DETERMINED', 'NONE_DETECTED', 'CONFLICT', 'NOT_ASSESSED',
                 'AGREES', 'SUSPECT', 'DISAGREES', 'UNAVAILABLE', 'NOT_REQUESTED', 'NOT_APPLICABLE',
                 'READS', 'PRIOR_DEPENDENT', 'READ_SUPPORTED', 'PRIOR_ONLY', 'UNRESOLVED',
                 'LABELS_AMBIGUOUS', 'GENOTYPE_UNRESOLVED', 'PROVISIONAL_LABEL'):
        assert f'`{term}`' in module, f'Missing assessment/outcome description: {term}'
    for filename in ('sequence_locus.tsv', 'sequence_blocks.tsv', 'sequence_haplotypes.tsv',
                     'sequence_haplotypes.fa', 'sample.diagnostics/', 'scripts/extract_genotype_haplotypes.py'):
        assert filename in module
    assert module.count('<!-- option-inventory:begin -->') == module.count('<!-- option-inventory:end -->') == 1
    recipe = re.search(r"printf '([^']+)' > limits.tsv", module)
    assert recipe and recipe.group(1).encode().decode('unicode_escape') == (
        'profile\tmy_limits\ntime_limit_seconds\t900\nfootprint_limit_bytes\t4294967296\n')
    assert 'never changes the prediction' in module and 'None is calibrated confidence.' in module
    algorithm = (root/'docs/algorithms/genotype.md').read_text()
    stages = re.findall(r'^### (\d+)\. ', algorithm, re.M)
    trace = re.findall(r'^(\d+)\. ', algorithm.split('## Worked trace')[1], re.M)
    assert stages == trace == [str(i) for i in range(1, 11)]
    print(f'Documentation: {len(pages)} pages, {links} local links; manual contract and ten-stage trace pass')


if __name__ == '__main__':
    check(pathlib.Path(sys.argv[1]).resolve())
