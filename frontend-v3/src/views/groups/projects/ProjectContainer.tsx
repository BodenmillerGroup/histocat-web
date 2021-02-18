import { Route, Switch, useParams, useRouteMatch } from "react-router-dom";
import React, { useEffect, useState } from "react";
import { useProjectsStore } from "modules/projects";
import { LoadingView } from "components/LoadingView";
import shallow from "zustand/shallow";
import { ProjectView } from "./ProjectView";
import { WebSocketManager } from "utils/WebSocketManager";

export function ProjectContainer() {
  const { path } = useRouteMatch();
  const { projectId } = useParams<{ projectId: string }>();
  const { setActiveProjectId, getProjectData } = useProjectsStore(
    (state) => ({
      setActiveProjectId: state.setActiveProjectId,
      getProjectData: state.getProjectData,
    }),
    shallow
  );
  const [loading, setLoading] = useState<boolean>(true);

  // Mounted
  useEffect(() => {
    const activeProjectId = +projectId;
    setActiveProjectId(activeProjectId);
    getProjectData(activeProjectId).then(() => setLoading(false));
    WebSocketManager.connect(activeProjectId);
  }, [setActiveProjectId, getProjectData]);

  // Unmounted
  useEffect(() => {
    return () => {
      setActiveProjectId(null);
      WebSocketManager.close();
    };
  }, []);

  if (loading) {
    return <LoadingView />;
  }

  return (
    <Switch>
      <Route exact={true} path={path} component={ProjectView} />
    </Switch>
  );
}
