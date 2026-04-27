workflow {
    // testing channel factory watchpath
    ch = channel.watchPath("${params.folderpath}*.pod5")
    ch.view { pod -> "File created or modified: $pod" }
}