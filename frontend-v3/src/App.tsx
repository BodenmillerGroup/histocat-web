import React, { Suspense, useEffect } from "react";
import { MainNavBar } from "views/MainNavBar";
import { FocusStyleManager } from "@blueprintjs/core";
import { LoginForm } from "views/auth/LoginForm";
import { BrowserRouter, Switch, Route } from "react-router-dom";
import { useAuthStore } from "modules/auth";
import history from "utils/history";
import { SignupForm } from "./views/auth/SignupForm";
import { PasswordRecoveryForm } from "./views/auth/PasswordRecoveryForm";
import { ResetPasswordForm } from "./views/auth/ResetPasswordForm";
import "./App.scss";
import { LoadingView } from "./components/LoadingView";
const GroupsView = React.lazy(() => import("./views/groups/GroupsView"));
const ProfileView = React.lazy(() => import("./views/profile/ProfileView"));
const AdminView = React.lazy(() => import("./views/admin/AdminView"));

FocusStyleManager.onlyShowFocusOnTabs();

function App() {
  const { loggedIn, checkLoggedIn } = useAuthStore();

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
    <BrowserRouter>
      <div className="bp3-dark">
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
    </BrowserRouter>
  );
}

export default App;
