#!/usr/bin/env python3

# Used for encoding image files to base64 for HTML embedding
from base64 import b64encode


def to_number(var):
    """Convert a string to an int or float if possible, otherwise return None.
       Used to auto-cast values from strings during QC parsing.
    """
    try:
        if '.' in var:
            raise ValueError
        return int(var)
    except ValueError:
        try:
            return float(var)
        except ValueError:
            return None


class QCLog(object):
    """
    Parses a QC log file (typically TSV format) into a Python dictionary.

    If no parser is provided, it assumes the file is tab-delimited and
    uses the first column as keys and the rest as values.
    """

    def __init__(self, log_file, parser=None):
        """
        Args:
            log_file: path to the log file or a pre-parsed dictionary
            parser: optional parser function for non-standard formats
        """
        self._log_file = log_file
        self._parser = parser
        self._dict = None
        self.__parse()

    def to_dict(self):
        """Return the parsed dictionary."""
        return self._dict

    def __parse(self):
        """Parse the log file (or use the parser if provided)."""
        if isinstance(self._log_file, dict):
            self._dict = self._log_file
            return

        if self._parser is None:
            # Default parser assumes TSV file
            d = {}
            with open(self._log_file, 'r') as fp:
                lines = fp.read().strip('\n')

            def to_number_or_str(var):
                """Try converting to number; fallback to string."""
                try:
                    if '.' in var:
                        raise ValueError
                    return int(var)
                except ValueError:
                    try:
                        return float(var)
                    except ValueError:
                        return var

            # Read line by line and populate the dict
            for i, line in enumerate(lines.split('\n')):
                arr = line.split('\t')
                if len(arr) == 1:
                    key = 'value' + str(i)
                    val = to_number(arr[0])
                elif len(arr) == 2:
                    key = arr[0]
                    val = to_number(arr[1])
                elif len(arr) > 2:
                    key = arr[0]
                    val = [to_number(v) for v in arr[1:]]
                else:
                    continue
                d[key] = val

            self._dict = d
        else:
            # Use a custom parser if provided
            self._dict = self._parser(self._log_file)


class QCPlot(object):
    """
    Encodes an image file to base64 and formats it as an embeddable
    HTML <img> element inside a <figure> tag.
    """

    def __init__(self, plot_file, caption=None, size_pct=100):
        self._plot_file = plot_file
        self._caption = caption
        self._size_pct = size_pct
        self._encoded = None
        self._img_type = None
        self.__encode()

    def to_html(self):
        """Return HTML figure with image and caption."""
        html = '''
        <figure style="display:inline-block">
          <img src="data:image/{img_type};base64,{encoded}" alt="{caption}" height="{size_pct}%"/>
          <figcaption style="text-align:center">{caption}</figcaption>
        </figure>
        '''
        return html.format(
            img_type=self._img_type,
            size_pct=self._size_pct,
            encoded=self._encoded,
            caption='' if self._caption is None else self._caption)

    def __encode(self):
        """Read and base64-encode the plot file."""
        self._encoded = b64encode(
            open(self._plot_file, 'rb').read()).decode("utf-8")
        self._img_type = os.path.splitext(self._plot_file)[-1].lstrip('.')


class QCCategory(object):
    """
    Organizes QC logs and plots into categories, supports nesting of subcategories.
    This allows structured and hierarchical QC reporting.
    """

    def __init__(self, cat_name, html_head='', html_foot='',
                 parser=None, map_key_desc=None, parent=None):
        """
        Args:
            cat_name: Category name for this QC section
            html_head: Optional HTML header (e.g., <h2> tag)
            html_foot: Optional HTML footer (e.g., line break)
            parser: Optional default parser for child QC logs
            map_key_desc: Dictionary mapping log keys to descriptions
            parent: Optional parent QCCategory (to create hierarchy)
        """
        self._cat_name = cat_name
        self._html_head = html_head
        self._html_foot = html_foot
        self._parser = parser
        self._map_key_desc = map_key_desc
        self._qc_logs = {}          # Keyed log files
        self._qc_plots = {}         # Keyed image plots
        self._child_categories = [] # Nested child categories
        if parent is not None:
            parent.add_category(self)

    def add_category(self, qc_category):
        """Add a nested QCCategory as a child."""
        self._child_categories.append(qc_category)

    def add_log(self, log_file, key=None):
        """Add a QCLog entry."""
        assert(key not in self._qc_logs)
        self._qc_logs[key] = QCLog(log_file, parser=self._parser)

    def add_plot(self, plot_file, key=None, caption=None, size_pct=100):
        """Add a QCPlot entry."""
        assert(key not in self._qc_plots)
        self._qc_plots[key] = QCPlot(
            plot_file,
            caption=caption if caption is not None else key,
            size_pct=size_pct
        )

    def to_dict(self):
        """Convert all logs and subcategories into a nested dictionary."""
        d = {}
        for key, qc_log in self._qc_logs.items():
            d_ = qc_log.to_dict()
            if len(d_) > 0:
                d[key] = d_
        for cat in self._child_categories:
            d_ = cat.to_dict()
            if len(d_) > 0:
                d[cat._cat_name] = d_
        return d

    def to_html(self):
        """Render the category and its contents into HTML."""
        html = ''
        html += self.__qc_logs_to_html()
        html += self.__qc_plots_to_html()
        for cat in self._child_categories:
            html += cat.to_html()
        if html == '':
            return ''
        else:
            return self._html_head + html + self._html_foot

    def __qc_logs_to_html(self):
        """Render QC logs into an HTML table."""
        if len(self._qc_logs) == 0:
            return ''

        html = '<table border="1" style="border-collapse:collapse">'

        # Build table header row
        header = '<tr><th bgcolor="#EEEEEE">'
        arr = [' ']
        if len(self._qc_logs) == 1 and list(self._qc_logs.keys())[0] is None:
            arr += ['Description']
        else:
            arr += self._qc_logs.keys()
        header += '</th><th bgcolor="#EEEEEE">'.join(arr) + '</th></tr>\n'

        # Collect all possible keys from all logs
        all_keys = {}
        for qc_log_k, qc_log_val in self._qc_logs.items():
            new_keys = dict.fromkeys(
                k for k in qc_log_val.to_dict() if k not in all_keys.keys())
            for new_key in new_keys:
                all_keys[new_key] = None

        # Build content rows
        content = ''
        for key in all_keys.keys():
            # Use description mapping if provided
            long_key_name = self._map_key_desc.get(key, key) if self._map_key_desc else key
            content += '<tr><th bgcolor="#EEEEEE" style="text-align:left">'\
                '{}</th><td>'.format(long_key_name)

            # Fill in row values across all logs
            qc_log_content = []
            for qc_log_key, qc_log in self._qc_logs.items():
                d = qc_log.to_dict()
                val = str(d.get(key, 'N/A'))
                qc_log_content.append(val)

            content += '</td><td>'.join(qc_log_content)
            content += '</td></tr>\n'

        html += header
        html += content
        html += '</table><br>\n'
        return html

    def __qc_plots_to_html(self):
        """Render QC plots as HTML figures."""
        html = ''
        for k, qc_plot in self._qc_plots.items():
            html += qc_plot.to_html()
        return html

