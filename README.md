# sequana-wrappers-lite

## :question: What is it 
a light version of the sequana-wrappers repository https://github.com/sequana/sequana-wrappers

## :dart: targetted audience

Sequana's developers.

The idea is that from sequana_pipetools, we can fetch this repo that will be a copy of the wrappers of sequana-wrappers without the test, documentation, large .git directory, etc

Yet, we need to make sure this is synchronised with the official wrappers, including tags.

Here we provide a simple utility update.py that fetches a tag and copy only the wrappers and environment.yaml files for you.

Only tagged version from sequana-wrappers should be included here.

From python using the update.py file provided:

    import update
    w = update.Wrapper()
    w.run("v24.8.29")

exit, and do:

    cd wrappers
    git commit  -m "Add wrapper.py files for tag v24.8.29"
    git tag v24.8.29
    git push origin main
    git push origin v24.8.29

Or even simpler, use

    python update.py --tag vX.Y.Z

and follow the instructions.

We already commited the older tag but for information and book-keeping, keep in mind the following issue. If you consider adding an older tag, this rewrite history of the original wrappers. This is not a big deal since this repo is not for developement. Yet, would be nice to have the same behaviour as the original wrappers if we do a git clone. So consider reset the HEAD to the commit of the newest tag whenever required using

    git checkout main
    git reset --hard v24.8.29 
    git push origin main --force


## 🗨️ Contacts <a name="contacts"></a>

Please fill an issue. 
