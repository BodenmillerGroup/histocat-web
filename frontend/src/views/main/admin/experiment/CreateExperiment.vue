<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="headline primary--text">Create Experiment</div>
      </v-card-title>
      <v-card-text>
        <template>
          <v-form v-model="valid" ref="form" lazy-validation>
            <v-text-field label="Name" v-model="name" :rules="nameRules"></v-text-field>
            <v-text-field label="Description" v-model="description"></v-text-field>
            <v-combobox
              v-model="tags"
              :filter="filter"
              :hide-no-data="!search"
              :items="items"
              :search-input.sync="search"
              hide-selected
              label="Tags"
              multiple
              small-chips
              solo
            >
              <template v-slot:no-data>
                <v-list-item>
                  <span class="subtitle-1">Create</span>
                  <v-chip label small>
                    {{ search }}
                  </v-chip>
                </v-list-item>
              </template>
              <template v-slot:selection="{ attrs, item, parent, selected }">
                <v-chip v-if="item === Object(item)" v-bind="attrs" :input-value="selected" label small>
                  <span class="pr-2">
                    {{ item.text }}
                  </span>
                  <v-icon small @click="parent.selectItem(item)">mdi-close </v-icon>
                </v-chip>
              </template>
              <template v-slot:item="{ index, item }">
                <v-text-field
                  v-if="editing === item"
                  v-model="editing.text"
                  autofocus
                  flat
                  background-color="transparent"
                  hide-details
                  solo
                  @keyup.enter="edit(index, item)"
                ></v-text-field>
                <v-chip v-else dark label small>
                  {{ item.text }}
                </v-chip>
                <v-spacer></v-spacer>
                <v-list-item-action @click.stop>
                  <v-btn icon @click.stop.prevent="edit(index, item)">
                    <v-icon>{{ editing !== item ? "mdi-pencil" : "mdi-check" }}</v-icon>
                  </v-btn>
                </v-list-item-action>
              </template>
            </v-combobox>
          </v-form>
        </template>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn @click="submit" :disabled="!valid">
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { IExperimentCreate } from "@/modules/experiment/models";
import { required } from "@/utils/validators";
import { Component, Vue, Watch } from "vue-property-decorator";

@Component
export default class CreateExperiment extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);

  readonly nameRules = [required];

  valid = false;
  name: string = "";
  description: string = "";

  editing = null;
  index = -1;
  nonce = 1;
  tags: any[] = [];
  search = null;

  get items(): any[] {
    const list = this.experimentContext.getters.tags;
    return list.map(item => {
      return {
        text: item
      };
    });
  }

  async mounted() {
    await this.experimentContext.actions.getTags();
  }

  reset() {
    this.name = "";
    this.description = "";
    this.tags = [];
    (this.$refs.form as any).resetValidation();
  }

  cancel() {
    this.$router.back();
  }

  @Watch("tags")
  onModelChange(val: any[], prev: any[]) {
    if (val.length === prev.length) {
      return;
    }

    this.tags = val.map(v => {
      if (typeof v === "string") {
        v = {
          text: v
        };
        this.items.push(v);
        this.nonce++;
      }
      return v;
    });
  }

  edit(index: number, item) {
    if (!this.editing) {
      this.editing = item;
      this.index = index;
    } else {
      this.editing = null;
      this.index = -1;
    }
  }

  filter(item, queryText: string, itemText: string) {
    const hasValue = val => (val != null ? val : "");

    const text = hasValue(itemText);
    const query = hasValue(queryText);

    return (
      text
        .toString()
        .toLowerCase()
        .indexOf(query.toString().toLowerCase()) > -1
    );
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const params: IExperimentCreate = {
        name: this.name
      };
      if (this.description) {
        params.description = this.description;
      }
      if (this.tags.length > 0) {
        params.tags = this.tags.map(tag => tag.text);
      }
      await this.experimentContext.actions.createExperiment(params);
      this.$router.back();
    }
  }
}
</script>
