# SYSADMIN.md: operational runbook for the ssign hosted service

Runbook for the lab member maintaining the ssign hosted web service after the
original author (Teo Reid) leaves Imperial. Written for a moderately
technical lab member, not a professional sysadmin.

> **Scope:** public web service only. The ssign pipeline (what you
> `pip install` or `docker run`) is documented under [`docs/`](docs/).
>
> **Status (2026-04-22):** skeleton. Sections marked "to be filled in" will
> be expanded as the hosted service is built out (post-publication,
> separate plan).

---

## 0. Who owns what, who to contact

- **Server owner / lab contact:** Dr. Sonja Billerbeck, `s.billerbeck@imperial.ac.uk`
- **Named sysadmin in the lab** (designated before Teo departs): _to be
  assigned_
- **Imperial RCS** (for network / firewall / subdomain issues):
  `https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/`
- **Domain registrar** (if we own our own domain): _to be filled in_
- **Cloud services we depend on** (beyond Imperial):
  - Let's Encrypt (SSL); free, automated
  - GitHub (source)
  - Zenodo (model weights + database deposits)
  - Docker Hub / GHCR (container image)

---

## 1. One-page "something is broken" checklist

Before digging in, check these in order. Most outages fall here.

1. **Is the machine on?** Power cycle if unreachable.
2. **Can you ping the server?** `ping ssign.<domain>`
3. **Is the service running?** SSH in and run `systemctl status ssign`.
4. **Is Docker up?** `docker ps` should list the ssign container.
5. **Is the SSL cert current?** `certbot certificates` (should show > 0 days
   remaining).
6. **Is the disk full?** `df -h`; if /var is >95% it can break everything.
7. **Are there OOM events?** `dmesg | grep -i "killed process"`.
8. **Is there a GitHub / Docker Hub outage?** Check their status pages before
   blaming yourself.

If none of the above, see §6 "When it really is broken."

---

## 2. Hardware + OS

- **Target hardware:** Lenovo ThinkStation P520 or Dell Precision T5820 (see
  resource-plan memo).
- **OS:** Ubuntu LTS (latest supported at deploy time; pin exact release).
- **Automatic security updates:** `unattended-upgrades` enabled for security
  patches only; auto-reboot at 03:00 local if kernel update demands it.
  Service comes back via `systemctl enable ssign`.
- **Full package upgrades:** _do not_ run `apt upgrade` / `do-release-upgrade`
  on a schedule. Do it intentionally, in a maintenance window, after taking a
  backup.
- **Backups:** `/etc`, `/var/lib/ssign`, user-submitted inputs and outputs
  under `/srv/ssign/submissions/`. Target backup location to be filled in.

_Config files + exact command list to be added when hardware is deployed._

---

## 3. SSL / HTTPS

- **Certificate authority:** Let's Encrypt (nonprofit, free).
- **Client:** `certbot` (standard Ubuntu package).
- **Renewal:** automatic via systemd timer, runs twice daily, renews when
  ≤30 days to expiry.
- **Email alerts on renewal failure:** yes, via `certbot register --email <addr>`.
  Failure emails should go to at least two lab members.
- **Cert lifetime:** 90 days. Renewal should happen silently ~every 60 days.
- **If renewal fails for >2 weeks:** site will show "Not Secure" warnings to
  users. Diagnose with `certbot renew --dry-run`; most causes are port 80
  blocked or DNS misconfigured.
- **Imperial subdomain option:** if we move to an `*.imperial.ac.uk`
  subdomain, Imperial ICT handles certs and this section can mostly be
  deleted.

_Exact config files + troubleshooting commands to be added once SSL is set up._

---

## 4. Networking

- **Ports open:** 80 (HTTP → redirects to HTTPS), 443 (HTTPS). Nothing else
  should be exposed to the public internet.
- **SSH:** behind Imperial VPN or a private tunnel (Tailscale / Imperial
  SSH bastion); never exposed to the public internet directly.
- **Firewall:** `ufw` or `nftables`, default-deny for inbound.

_Specific firewall rules + VPN/bastion setup to be added at deployment time._

---

## 5. ssign service

- **Runs as:** `systemd` unit named `ssign.service`, running a Docker
  container (`docker run billerbeck-lab/ssign:<version>`).
- **Log location:** `journalctl -u ssign` (systemd) and
  `/var/log/ssign/` (application).
- **Data volumes:** databases at `/srv/ssign/databases` (read-only), user
  submissions at `/srv/ssign/submissions` (writable, auto-pruned after 30
  days).
- **Upgrade procedure:** pull new image, `systemctl restart ssign`. Roll back
  by pinning the previous image tag and restarting.

_Exact service file, Docker-compose, and data-retention cron to be added at
deployment time._

---

## 6. When it really is broken

1. **Gather info first.** Don't `reboot` as a reflex.
   - `systemctl status ssign` (service state)
   - `journalctl -u ssign --since "1 hour ago" | tail -100` (recent errors)
   - `docker logs ssign` (container output)
   - `df -h` / `free -h` (disk / memory pressure)
   - `dmesg | tail -50` (kernel warnings)
2. **Take the service offline** if it's misbehaving visibly.
   `systemctl stop ssign`. Users see "service unavailable" during investigation.
3. **Roll back** if the last change was a new Docker image:
   `docker pull billerbeck-lab/ssign:<previous-tag>` and restart.
4. **Restore from backup** if data is corrupted.
5. **Ask for help.** Contact list in §0. If Imperial RCS is involved,
   reference the original server-provisioning ticket.
6. **Document** what broke, what you did, and what worked. Append to this
   file (or the repo's issue tracker) so the next person benefits.

---

## 7. Known long-term risks

From the longevity memo (`project_longevity_commitment.md`):

- **Hardware MTBF** ~5 years for refurbished workstations. Plan hardware
  replacement in year 4–5.
- **DTU licenses** (SignalP 6.0, DeepLocPro) for in-container redistribution:
  policy may change. Fallback is "user installs separately," already
  documented in the install docs.
- **External service changes** (Docker Hub pulls, Zenodo URLs): we mirror
  critical artifacts to institutional RDS as backup.
- **Upstream tool updates** may invalidate frozen results. The frozen v1.0.0
  image exists precisely so this doesn't force you to update.

---

## 8. If you want to decommission the service

Acceptable. The publication + GitHub repo + Docker image + Zenodo deposits
stand alone; anyone can still run ssign locally. Steps:

1. Post a notice on the webservice landing page ("hosted service
   decommissioned YYYY-MM-DD; please run locally via GitHub / Docker").
2. Leave the notice up for ≥90 days before taking the site down.
3. Update the README in the GitHub repo to point users at local install
   instructions.
4. Power down the hardware. Keep the databases + Docker image archived on
   institutional RDS for 2+ years.
5. Email Sonja + Imperial RCS so they know the service is gone.

The paper remains valid. Users lose convenience; they do not lose access
to the tool.

---

_This runbook will be expanded as the hosted service is deployed. Last updated:
2026-04-22 (skeleton only)._
