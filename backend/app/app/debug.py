import ptvsd


def init_debug():
    # Allow other computers to attach to ptvsd at this IP address and port.
    ptvsd.enable_attach(address=('172.26.0.8', 5678), redirect_output=True)

    # Pause the program until a remote debugger is attached
    # ptvsd.wait_for_attach()
