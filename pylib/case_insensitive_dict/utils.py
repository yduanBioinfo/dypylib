import unicodedata


# Encode/Decode
def clean_unicode(value):
    try:
        return unicodedata.normalize(
            'NFKD', value).encode('ascii', 'ignore').decode('ascii')
    except Exception:
        return value
