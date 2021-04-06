import React, { Suspense, useEffect } from "react";
import { MainNavBar } from "views/MainNavBar";
import { Classes, FocusStyleManager } from "@blueprintjs/core";
import { LoginForm } from "views/auth/LoginForm";
import { Router, Switch, Route } from "react-router-dom";
import { useAuthStore } from "modules/auth";
import history from "utils/history";
import { SignupForm } from "./views/auth/SignupForm";
import { PasswordRecoveryForm } from "./views/auth/PasswordRecoveryForm";
import { ResetPasswordForm } from "./views/auth/ResetPasswordForm";
import { LoadingView } from "./components/LoadingView";
import shallow from "zustand/shallow";

import "react-grid-layout/css/styles.css";
import "react-resizable/css/styles.css";
// import "vitessce/css/index.scss";
import "./App.scss";
const GroupsView = React.lazy(() => import("./views/groups/GroupsView"));
const ProfileView = React.lazy(() => import("./views/profile/ProfileView"));
const AdminView = React.lazy(() => import("./views/admin/AdminView"));

FocusStyleManager.onlyShowFocusOnTabs();

function App() {
  const { loggedIn, checkLoggedIn } = useAuthStore(
    (state) => ({
      loggedIn: state.loggedIn,
      checkLoggedIn: state.checkLoggedIn,
    }),
    shallow
  );

  useEffect(() => {
    checkLoggedIn();
  }, [checkLoggedIn]);

  if (!loggedIn) {
    switch (history.location.pathname) {
      case "/signup":
        return <SignupForm />;
      case "/password-recovery":
        return <PasswordRecoveryForm />;
      case "/reset-password":
        return <ResetPasswordForm />;
      default:
        return <LoginForm />;
    }
  }

  return (
    <Router history={history}>
      <div className={Classes.DARK}>
        <MainNavBar />
        <main>
          <Switch>
            <Route path="/profile">
              <Suspense fallback={<LoadingView />}>
                <ProfileView />
              </Suspense>
            </Route>

            <Route path="/admin">
              <Suspense fallback={<LoadingView />}>
                <AdminView />
              </Suspense>
            </Route>

            <Route path="/">
              <Suspense fallback={<LoadingView />}>
                <GroupsView />
              </Suspense>
            </Route>
          </Switch>
        </main>
      </div>
    </Router>
  );
}

export default App;
