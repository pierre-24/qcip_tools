from collections import ChainMap


class DispatcherMeta(type):
    def __new__(mcs, name, bases, attrs):
        callbacks = ChainMap()

        # extend parent callbacks
        maps = callbacks.maps
        for base in bases:
            if isinstance(base, DispatcherMeta):
                maps.extend(base.__callbacks__.maps)

        # add callbacks dictionary and dispatcher
        attrs['__callbacks__'] = callbacks
        attrs['dispatcher'] = property(lambda obj: callbacks)
        cls = super().__new__(mcs, name, bases, attrs)
        return cls

    def set_callback(cls, key, callback):
        """register callback

        :param key: key
        :type key: str
        :param callback: callback function
        :type callback: function
        """

        cls.__callbacks__[key] = callback
        return callback

    def register(cls, key):
        """decorator @Class.register(key)

        :param key: key
        :type key: str
        """
        def wrapper(callback):
            return cls.set_callback(key, callback)
        return wrapper


class Dispatcher(metaclass=DispatcherMeta):
    """Implement the dispatcher pattern:

    + use ``@Class.register(key)`` to register a new callback under  ``key`` ;
    + use ``dispatch(key)`` to access this callback.

    Source: https://zestedesavoir.com/tutoriels/1226/le-pattern-dispatcher-en-python/ (in french).
    """

    def dispatch(self, key, default=None):
        """Access to a registered callback

        :param key: the key
        :type key: str
        :param default: default value
        :rtype: function
        """

        return self.dispatcher.get(key, default)
