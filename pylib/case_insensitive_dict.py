# Encode/Decode
def clean_unicode(value):
    try:
        return unicodedata.normalize(
            'NFKD', value).encode('ascii', 'ignore').decode('ascii')
    except Exception:
        return value

class CaseInsensitiveDict(dict):
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        if isinstance(key, str):
            key = clean_unicode(key.lower())
        if isinstance(value, dict):
            value = CaseInsensitiveDict(value)
        if isinstance(value, list):
            value = [
                CaseInsensitiveDict(k)
                if isinstance(k, dict) else k for k in value]
        super(CaseInsensitiveDict, self).__setitem__(key, value)

    def __getitem__(self, key):
        if isinstance(key, str):
            return super(CaseInsensitiveDict, self).__getitem__(key.lower())
        return super(CaseInsensitiveDict, self).__getitem__(key)

    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except Exception:
            return default

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError("Update expected at most 1 arguments, got {}"
                                .format(len(args)))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]
