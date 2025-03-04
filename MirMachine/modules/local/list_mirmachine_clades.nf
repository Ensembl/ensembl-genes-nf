process LIST_MIRMACHINE_CLADES {
    label 'mirMachine'
    
    output:
    path "clades.txt", emit: clades

    script:
    """
    MirMachine.py --print-all-nodes > clades.txt
    """
}