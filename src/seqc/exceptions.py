class RetryLimitExceeded(Exception):
    pass


class InstanceNotRunningError(Exception):
    pass


class EC2RuntimeError(Exception):
    pass


class ConfigurationError(Exception):
    pass


class ArgumentParserError(Exception):
    pass


class EmptyMatrixError(Exception):
    pass
