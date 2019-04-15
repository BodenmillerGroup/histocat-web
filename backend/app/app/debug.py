import ptvsd


def init_debug():
    # Allow other computers to attach to ptvsd at this IP address and port.
    # ptvsd.enable_attach(address=('192.168.0.9', 5678), redirect_output=True)
    ptvsd.enable_attach(redirect_output=True)

    # Pause the program until a remote debugger is attached
    # ptvsd.wait_for_attach()
