const env = process.env.FRONTEND_ENV;

let envApiUrl = "";
let appNameTmp = process.env.VUE_APP_NAME;

if (env === "production") {
  envApiUrl = `http://${process.env.VUE_APP_DOMAIN_PROD}`;
} else if (env === "staging") {
  envApiUrl = `http://${process.env.VUE_APP_DOMAIN_STAG}`;
  appNameTmp = `${appNameTmp} (test)`;
} else {
  envApiUrl = `http://${process.env.VUE_APP_DOMAIN_DEV}`;
  appNameTmp = `${appNameTmp} (dev)`;
}

export const apiUrl = `${envApiUrl}/api/v1`;
export const appName = appNameTmp;
