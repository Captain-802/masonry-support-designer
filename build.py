"""Inlines the engines and the data into the page template (one HTML file, no external resources).

In the repository (layout src/template.html, src/sms.js, src/ec3_adapter.js, src/beam-v03-engine.js, src/v03.css,
src/data.json) this file is build.py:
    python build.py        ->  index.html
"""
import os


def _js(text):
    return text.replace('</script', '<\\/script')


def make_pages(template, engine, data_json, beamv03='', adapter='', v03css=''):
    """-> (body, full): the page without and with the <html> skeleton."""
    assert '/*@@SMS_JS@@*/' in template and '@@DATA@@' in template
    body = (template.replace('/*@@SMS_JS@@*/', _js(engine))
            .replace('/*@@BEAMV03_JS@@*/', _js(beamv03))
            .replace('/*@@EC3_ADAPTER_JS@@*/', _js(adapter))
            .replace('/*@@V03_CSS@@*/', v03css)
            .replace('@@DATA@@', data_json.replace('</', '<\\/')))
    full = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n' + body + '\n</html>\n')
    return body, full


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))

    def read(*p):
        f = os.path.join(here, *p)
        return open(f, encoding='utf-8').read() if os.path.exists(f) else ''

    _, full = make_pages(read('src', 'template.html'), read('src', 'sms.js'), read('src', 'data.json'),
                         read('src', 'beam-v03-engine.js'), read('src', 'ec3_adapter.js'), read('src', 'v03.css'))
    open(os.path.join(here, 'index.html'), 'w', encoding='utf-8', newline='\n').write(full)
    print('index.html', len(full) // 1024, 'KB')
