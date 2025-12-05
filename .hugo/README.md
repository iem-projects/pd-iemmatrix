iemmatrix web documentation
===========================

https://pd.iem.sh/iemmatrix/

# Workflow

This page uses [`hugo`](https://gohugo.io) to turn a bunch of markdown files (`.md`) into a static webpage.

1. install Hugo: https://gohugo.io/installation/
2. Clone this repository
3. Open a terminal, navigate into the cloned repository
4. change into the (hidden) `.hugo` directory
   ```sh
   cd .hugo
   ```
4. Run `hugo server` to startup a live renderer.
   You should see something like this:

5. In your browser, open http://localhost:1313/ to get a live preview
   (the actual port number (here `1313` may be different, depending on port availability on your system); it is shown at the end of the output when starting `hugo server`)

6. Edit the `.md` files in the `documentation/` directory (*outside* of the `.hugo` directory).
   As soon as you <kbd>save</<kbd> a Markdown file, the live preview will be updated.
7. Once finished, stop the live renderer with <kbd>Ctrl</kbd>-<kbd>C</kbd> (as shown in the startup output)


# rebuilding the public webpage

To rebuild the public webpage of iemmatrix (available at https://software.iem.at/iemmatrix)
just push your changes to the canonical repository (https://git.iem.at/pd/iemmatrix).


## *only* rebuild the public webpage

When pushing *anything* to the canonical repository,
a number of CI jobs are run (e.g. building/signing/notarizing `iemmatrix` for macOS, Linux & Windows,... ).

If you know that your changeset only contains changes to the HTML documentation,
you might want to skip the re-compilation jobs.

To only run the `pages` job (which creates the public HTML documentation), push your changes like this:

```sh
git push -o ci.variable="IEM_CI_JOBS=pages"
```


Alternatively, if you want to skip the CI entirely, you can use:
```sh
git push -o ci.skip
```
