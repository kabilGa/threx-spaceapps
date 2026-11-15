# Submit checklist

The local git commit is already made. Nothing has been pushed.

## A. Push to GitHub (if the GitHub CLI `gh` is installed and `gh auth status` says you are logged in)
```bash
gh repo create kabilGa/threx-spaceapps --public --source=. --push
```

## B. Otherwise (no `gh`)
1. On https://github.com/new create an EMPTY public repository named `threx-spaceapps` (no README, no licence, no .gitignore).
2. In the `threx-spaceapps` folder:
```bash
git push -u origin main
```
(The `origin` remote is already set. If Git asks you to log in, use the browser sign-in.)

Then open https://github.com/kabilGa/threx-spaceapps/blob/main/slides/threx_deck.pdf and check the deck displays.

## C. Fill in the Space Apps form
Open the challenge page and your project form. Copy `form_text.md` into the form field by field, in its order: Project name, High-level summary, Demo link, Link to final project, Detailed description, NASA data used, Other data and tools, References, Use of Artificial Intelligence. Upload the four images listed at the end of `form_text.md`, tick the two boxes, then click **Submit for Judging** before Sunday 15 November 2026, 23:59 Tunis time.

## D. Optional: public web app
On https://share.streamlit.io: repository `kabilGa/threx-spaceapps`, branch `main`, main file `app.py`, Python 3.14 (the version the pinned `requirements.txt` was tested with; 3.13 may also work). The submission does not need it: the app also runs locally with `streamlit run app.py`.
