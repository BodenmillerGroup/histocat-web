import React from "react";
import TitleInfo from "vitessce/components/TitleInfo";
import { SlidesView } from "./SlidesView";
import { CoordinationScopes } from "vitessce/types";

type SlidesViewSubscriberProps = {
  coordinationScopes: CoordinationScopes;
  removeGridComponent(event: any): void;
};

export function SlidesViewSubscriber(props: SlidesViewSubscriberProps) {
  const { removeGridComponent } = props;
  return (
    <TitleInfo
      title="Slides"
      removeGridComponent={removeGridComponent}
      isScroll={true}
      theme="dark"
      isReady={true}
    >
      <SlidesView />
    </TitleInfo>
  );
}
