"""Functions to calculate features associated with each amino acid in a protein sequence."""
import re

# Helper functions
def get_window(seq, pos, window_size):
    """Return window of length window_size centered at position pos in seq.

    If the window exceeds the bounds of the seq, get_window returns
    the maximal possible window. Thus, the window at the upper and
    lower bounds are actually right- and left- facing half windows,
    respectively.

    Parameters
    ----------
        seq : string
            Protein sequence as string.
        pos : int
            Index of the center position of the window.
        window_size : int
            Total number of symbols in window, including the center
            symbol. Must be an odd number.

    Returns
    -------
        window : string
            Window of length window_size centered at position pos in
            seq.
    """
    if pos < 0 or pos > len(seq) - 1:
        raise ValueError('Pos is outside the bounds of seq.')
    if window_size % 2 == 0:
        raise ValueError('Size is even.')
    delta = (window_size - 1) // 2

    lower = pos - delta
    if lower < 0:
        lower = 0
    upper = pos + delta + 1  # Add 1 to include upper bound
    if upper > len(seq):
        upper = len(seq)
    return seq[lower:upper]


def get_X_count(seq, X):
    """Return count of symbols in X in seq.

    Parameters
    ----------
        seq : string
            Protein sequence as string.
        X : string or list
            Symbols to count as string or list.

    Returns
    -------
        X_count : int
            Count of symbols in X in seq.
    """
    X_count = 0
    for sym in seq:
        if sym in X:
            X_count += 1
    return X_count


def get_X_fractions(seq, X, window_size):
    """Return fractions of symbols in X in sliding window across seq.

    Parameters
    ----------
        seq : string
            Protein sequence as string.
        X : string or list
            Symbols to count as string or list.
        window_size : int
            Total number of symbols in window, including the center
            symbol. Must be an odd number.

    Returns
    -------
        X_fractions : list
            Fractions of symbols in X in sliding window across seq
    """
    X_fractions = []
    for i in range(len(seq)):
        window = get_window(seq, i, window_size)
        X_count = get_X_count(window, X)
        X_fractions.append(X_count / len(window))
    return X_fractions


def get_regex_count(seq, regex):
    """Return count of patterns matching regex in seq."""
    return len(re.findall(regex, seq))


def get_regex_fractions(seq, regex, window):
    """Return fractions of patterns matching X in sliding window across seq."""
    pass


def get_hydrophobicity(seq):
    """
    Return average hydrophobicity of symbols in X in seq.

    Parameters
    ----------
        seq : string
            Protein sequence as string.

    Returns
    -------
        hydrophobicity : int
            Score of hydrophobicity with most hydrophobic at 1, and most
            hydrophilic at 0.
    """
    hydrophobicity_dict = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
                           'M': 1.9, 'A': 1.8, 'W': -0.9, 'G': -0.4, 'T': -0.7,
                           'S': -0.8, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'N': -3.5,
                           'D': -3.5, 'Q': -3.5, 'E': -3.5, 'K': -3.9, 'R': -4.5}
    for key in hydrophobicity_dict.keys():
        hydrophobicity_dict[key] = (hydrophobicity_dict[key] + 4.5) / 9
    seq_hydrophobicities = [hydrophobicity_dict.get(N, 0) for N in seq]
    return sum(seq_hydrophobicities) / len(seq_hydrophobicities)


def get_polarity(seq):
    """
    Return average polarity of symbols in X in seq.

    Parameters
    ----------
        seq : string
            Protein sequence as string.

    Returns
    -------
        polarity : int
            Score of average polarity ranging from 0 to 1, where 1 is all
            polar and 0 is all nonpolar.
    """
    polarity_dict = {'I': 0, 'V': 0, 'L': 0, 'F': 0, 'C': 0, 'M': 0, 'A': 0,
                     'W': 0, 'G': 0, 'T': 1, 'S': 1, 'Y': 1, 'P': 0, 'H': 1,
                     'N': 1, 'D': 1, 'Q': 1, 'E': 1, 'K': 1, 'R': 1}
    seq_polarities = [polarity_dict.get(N, 0.5) for N in seq]
    return sum(seq_polarities) / len(seq_polarities)


def get_content_scores(seq, content_func, window_size):
    """
    Return amino acid content score (scaled to 0-1) of amino acids in seq.

    Parameters
    ----------
        seq : string
            Protein sequence as string.
        content_func : function
            Gives content score for a specified sequence, ranging from 0 to 1.
        window_size : int
            Total number of symbols in window, including the center
            symbol. Must be an odd number.

    Returns
    -------
        content_scores : list
            Scores of content in sliding window across seq
    """
    content_scores = []
    for i in range(len(seq)):
        window = get_window(seq, i, window_size)
        score = content_func(seq)
        content_scores.append(score)
    return content_scores


# TODO: Not all features will necessarily accept a window_size parameter, so we may need to rework this
def get_features(seq, features, window_size):
    """Return values for each feature in features at each symbol in seq."""
    label2idx = {}  # Dictionary of (lists of feature values) keyed by feature label
    for feature_label, function in features.items():
        label2idx[feature_label] = function(seq, window_size)
    feature_labels = list(features)
    idx2label = []  # List of (dictionaries of feature values keyed by feature label)
    for i in range(len(seq)):
        idx2label.append({feature_label: label2idx[feature_label][i] for feature_label in feature_labels})
    return idx2label


# Feature functions and feature dictionary
features = {'hydrophobicity_content': lambda seq, window_size: get_content_scores(seq, get_hydrophobicity, window_size),
            'polarity_content': lambda seq, window_size: get_content_scores(seq, get_polarity, window_size),
            'aromatic_content': lambda seq, window_size: get_X_fractions(seq, 'FYWH', window_size)}
