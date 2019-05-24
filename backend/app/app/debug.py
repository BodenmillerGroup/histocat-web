def init_debug():
    # VS Code Debugging

    # Allow other computers to attach to ptvsd at this IP address and port.
    import ptvsd
    ptvsd.enable_attach(redirect_output=True)
    # Pause the program until a remote debugger is attached
    # ptvsd.wait_for_attach()


    # PyCharm Debugging
    # TODO: Don't forget to modify IP address!!

    # import pydevd_pycharm
    # pydevd_pycharm.settrace('130.60.106.83', port=5679, stdoutToServer=True, stderrToServer=True)
